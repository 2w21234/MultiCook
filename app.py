from flask import Flask, render_template, request, jsonify, Response, send_file
import os
import signal
import subprocess
import shutil
from werkzeug.utils import secure_filename
import re
import shlex
from subprocess import PIPE, Popen

app = Flask(__name__)

UPLOAD_FOLDER = 'uploads'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
process = None

@app.route('/')
def index():
    # Serve the main navigation page (base_navigation_template)
    return render_template('base_navigation_template.html')

@app.route('/imputation')
def imputation():
    # Serve the page for Step 1: Run Each Imputation
    return render_template('index.html')

@app.route('/merge_page')
def merge_page():
    # Serve the page for Step 2: Merge Results
    return render_template('merge.html')


"""
@app.route('/')
def index():
    return render_template('index.html')

@app.route('/merge_page')
def merge_page():
    return render_template('merge.html')
"""

@app.route('/merge_results', methods=['POST'])
def merge_results():
    try:
        working_directory = request.form.get('working_directory')
        if not working_directory:
            raise ValueError("Working directory not specified.")
        
        working_dir_path = os.path.join(UPLOAD_FOLDER, working_directory)
        output_path = os.path.join(working_dir_path, 'output')
        os.makedirs(output_path, exist_ok=True)

        input_list_path = os.path.join(output_path, 'input_list.txt')
        merge_command = [
            'python', 'Merge.py',
            '-i', input_list_path,
            '-o', os.path.join(output_path, 'Merge', 'result')
        ]
        
        print("Generated merge command:", ' '.join(merge_command))

        def generate():
            yield f"Executing command: {' '.join(merge_command)}\n\n"
            process = Popen(shlex.split(' '.join(merge_command)), stdout=PIPE, stderr=PIPE, text=True)
            for line in iter(process.stdout.readline, ''):
                yield line
            process.stdout.close()
            process.wait()

        return Response(generate(), mimetype='text/plain')

    except Exception as e:
        return jsonify(success=False, error=str(e))




@app.route('/upload_files', methods=['POST'])
def upload_files():
    try:
        working_directory = request.form.get('working_directory')
        if not working_directory:
            raise ValueError("Working directory not specified.")
        
        working_dir_path = os.path.join(UPLOAD_FOLDER, working_directory)
        output_path = os.path.join(working_dir_path, 'output')
        os.makedirs(output_path, exist_ok=True)

        # Get all file inputs and weights
        files = request.files.getlist('result_files[]')
        weight = request.form.get('weights[]')

        if not files or not weight:
            raise ValueError("Both files and weight must be provided.")

        # Get the next available additional folder number
        existing_folders = [d for d in os.listdir(output_path) 
                          if os.path.isdir(os.path.join(output_path, d)) 
                          and d.startswith('additional_')]
        next_number = 1 if not existing_folders else max([int(d.split('_')[1]) for d in existing_folders]) + 1

        # Create new additional folder
        additional_folder_name = f'additional_{next_number}'
        additional_folder_path = os.path.join(output_path, additional_folder_name)
        os.makedirs(additional_folder_path, exist_ok=True)

        # Process files
        if len(files) > 1:  # Multiple files selected
            # Save all files
            for file in files:
                if file.filename:
                    filename = secure_filename(file.filename)
                    file_save_path = os.path.join(additional_folder_path, filename)
                    file.save(file_save_path)
            
            # For input_list.txt, use only the first file after sorting
            filenames = [secure_filename(f.filename) for f in files if f.filename]
            sorted_filenames = sorted(filenames)
            if sorted_filenames:
                first_filename = sorted_filenames[0]
                first_file_path = os.path.join(additional_folder_path, first_filename)
                # Append to input_list.txt
                with open(os.path.join(output_path, 'input_list.txt'), 'a') as f:
                    f.write(f"{first_file_path} {weight}\n")
        else:  # Single file
            file = files[0]
            if file.filename:
                filename = secure_filename(file.filename)
                file_save_path = os.path.join(additional_folder_path, filename)
                file.save(file_save_path)
                # Append to input_list.txt
                with open(os.path.join(output_path, 'input_list.txt'), 'a') as f:
                    f.write(f"{file_save_path} {weight}\n")

        return jsonify(success=True)

    except Exception as e:
        print(f"Error in upload_files: {str(e)}")
        return jsonify(success=False, error=str(e))



@app.route('/run_imputation', methods=['POST'])
def run_imputation():
    try:
        global process
        print("Received POST request for imputation")

        # Handle input files
        input_files = request.files.getlist('input_files[]')
        if not input_files:
            raise ValueError("No input files uploaded.")

        # Retrieve working directory name
        working_directory = request.form.get('working_directory')
        if not working_directory:
            raise ValueError("Working directory not specified.")

        # Create working directory and input/output folders
        working_dir_path = os.path.join(UPLOAD_FOLDER, working_directory)
        input_folder = os.path.join(working_dir_path, 'input')
        output_folder = os.path.join(working_dir_path, 'output')
        os.makedirs(input_folder, exist_ok=True)
        os.makedirs(output_folder, exist_ok=True)

        input_file_prefix = None
        for file in input_files:
            filename = secure_filename(file.filename)
            if not input_file_prefix:
                input_file_prefix = os.path.splitext(filename)[0]  # Extract prefix from the first file
            file_path = os.path.join(input_folder, filename)
            file.save(file_path)
            print(f"Saved input file: {file_path}")

        if not input_file_prefix:
            raise ValueError("Could not determine input file prefix.")

        # Retrieve form data
        tools = request.form.getlist('tools[]')
        weights = request.form.getlist('weights[]')
        genome_assembly = request.form.get('genome_assembly', '')
        memory = request.form.get('memory', '')

        # Handle reference panel files
        ref_panels = request.files.getlist('ref_panels[]')
        if not ref_panels:
            raise ValueError("No reference panel files uploaded.")

        ref_panel_command_paths = []
        for tool in tools:
            if tool == "CookHLA":
                # Extract common prefix for CookHLA references
                cookhla_files = [file for file in ref_panels if any(ext in file.filename for ext in ['.bed', '.bim', '.fam', 'FRQ.frq', 'bgl.phased', 'markers'])]
                if len(cookhla_files) != 6:
                    raise ValueError("CookHLA requires exactly 6 reference panel files.")
                common_prefix = os.path.commonprefix([os.path.basename(file.filename) for file in cookhla_files])
                common_prefix = re.sub(r'\.[^.]*$', '', common_prefix)  # Remove extension if present
                ref_panel_command_paths.append(os.path.join(UPLOAD_FOLDER, 'reference', common_prefix))
            elif tool == "HIBAG":
                # Handle HIBAG reference file
                hibag_file = next((file for file in ref_panels if file.filename.endswith('.RData')), None)
                if not hibag_file:
                    raise ValueError("HIBAG requires a reference panel file with .RData extension.")
                hibag_path = os.path.join(UPLOAD_FOLDER, 'reference', secure_filename(hibag_file.filename))
                ref_panel_command_paths.append(hibag_path)

        # Save reference panel files
        for file in ref_panels:
            filename = secure_filename(file.filename)
            file_path = os.path.join(UPLOAD_FOLDER, 'reference', filename)
            os.makedirs(os.path.dirname(file_path), exist_ok=True)
            file.save(file_path)
            print(f"Saved reference panel file: {file_path}")

        # Generate command
        command = [
            'python', 'Imputation_single.py',
            '-i', os.path.join(input_folder, input_file_prefix),
            '-o', output_folder,
            '-t', ' '.join(tools),
            '-w', ' '.join(weights),
            '-hg', genome_assembly,
            '-m', memory,
            '-r', ' '.join(ref_panel_command_paths)
        ]

        print("Generated command:", ' '.join(command))

        # Generator function to execute command and stream output
        def generate():
            global process
            yield f"Executing command: {' '.join(command)}\n\n"
            process = Popen(shlex.split(' '.join(command)), stdout=PIPE, stderr=PIPE, text=True)
            for line in iter(process.stdout.readline, ''):
                yield line
            process.stdout.close()
            process.wait()

        return Response(generate(), mimetype='text/plain')

    except Exception as e:
        print(f"Exception occurred: {e}")
        return jsonify(error=str(e)), 500

@app.route('/download_folder/<path:folder_name>', methods=['GET'])
def download_folder(folder_name):
    try:
        # Ensure the folder path exists
        folder_path = os.path.join(UPLOAD_FOLDER, folder_name, 'output')
        if not os.path.exists(folder_path):
            return jsonify({"error": f"Folder '{folder_name}/output' does not exist."}), 404

        # Check if folder is empty
        folder_contents = os.listdir(folder_path)
        if not folder_contents:
            return jsonify({"error": f"Folder '{folder_name}/output' is empty."}), 400

        # Log folder contents for debugging
        print(f"Contents of '{folder_name}/output': {folder_contents}")

        # Create a zip file only if it doesn't already exist
        zip_file_path = f"{folder_path}.zip"
        if not os.path.exists(zip_file_path):
            shutil.make_archive(folder_path, 'zip', folder_path)

        # Serve the zip file
        return send_file(zip_file_path, as_attachment=True)

    except Exception as e:
        # Log the exception for debugging
        print(f"Error in /download_folder: {str(e)}")
        return jsonify({"error": "An error occurred while processing your request."}), 500

"""
@app.route('/run_imputation', methods=['POST'])
def run_imputation():
    try:
        # Get parameters from the request
        working_directory = request.form.get('working_directory')
        genome_assembly = request.form.get('genome_assembly')
        memory = request.form.get('memory')

        if not working_directory:
            raise ValueError("Working directory not specified")

        # Get lists of tools and weights
        tools = request.form.getlist('tools[]')
        weights = request.form.getlist('weights[]')
        
        # Get input files
        input_files = request.files.getlist('input_files[]')
        ref_panels = request.files.getlist('ref_panels[]')

        if not tools or not weights or not input_files or not ref_panels:
            raise ValueError("Missing required parameters")

        # Create working directory structure
        working_dir_path = os.path.join(UPLOAD_FOLDER, working_directory)
        os.makedirs(working_dir_path, exist_ok=True)

        # Save input files
        for file in input_files:
            if file.filename:
                filename = secure_filename(file.filename)
                file.save(os.path.join(working_dir_path, filename))

        # Process each tool
        for i, (tool, weight) in enumerate(zip(tools, weights)):
            # Create tool-specific directory
            tool_dir = os.path.join(working_dir_path, f'tool_{i+1}')
            os.makedirs(tool_dir, exist_ok=True)

            # Save reference panel files for this tool
            ref_panel = ref_panels[i] if i < len(ref_panels) else None
            if ref_panel and ref_panel.filename:
                ref_filename = secure_filename(ref_panel.filename)
                ref_panel.save(os.path.join(tool_dir, ref_filename))

        # Prepare imputation command
        command = [
            'python3', 'run_imputation.py',
            '--input', working_dir_path,
            '--genome-assembly', genome_assembly
        ]
        
        if memory:
            command.extend(['--memory', memory])

        print("Generated command:", ' '.join(command))

        def generate():
            yield f"Executing command: {' '.join(command)}\n\n"
            process = Popen(command, stdout=PIPE, stderr=PIPE, text=True)
            
            for line in iter(process.stdout.readline, ''):
                yield line
            
            process.stdout.close()
            process.wait()

        return Response(generate(), mimetype='text/plain')

    except Exception as e:
        return jsonify(success=False, error=str(e))


@app.route('/run_imputation', methods=['POST'])
def run_imputation():
    try:
        # Get parameters from the request
        working_directory = request.form.get('working_directory')
        marker = request.form.get('marker')
        model = request.form.get('model')
        window = request.form.get('window')

        if not all([working_directory, marker, model, window]):
            raise ValueError("Missing required parameters")

        # Construct paths
        working_dir_path = os.path.join(UPLOAD_FOLDER, working_directory)
        output_path = os.path.join(working_dir_path, 'output')

        # Prepare the command
        command = [
            'python3', 'CookHLA.py',
            '--input', os.path.join(working_dir_path, 'MHC'),
            '--hla-list', os.path.join(working_dir_path, 'HLA_ALLELE.txt'),
            '--model', model,
            '--window', window,
            '--out', os.path.join(output_path, 'CookHLA_OUT')
        ]

        if marker:
            command.extend(['--marker', marker])

        print("Generated command:", ' '.join(command))

        def generate():
            yield f"Executing command: {' '.join(command)}\n\n"
            process = Popen(command, stdout=PIPE, stderr=PIPE, text=True)
            
            for line in iter(process.stdout.readline, ''):
                yield line
            
            process.stdout.close()
            process.wait()

        return Response(generate(), mimetype='text/plain')

    except Exception as e:
        return jsonify(success=False, error=str(e))

@app.route('/upload_files', methods=['POST'])
def upload_files():
    try:
        working_directory = request.form.get('working_directory')
        if not working_directory:
            raise ValueError("Working directory not specified.")
        
        working_dir_path = os.path.join(UPLOAD_FOLDER, working_directory)
        output_path = os.path.join(working_dir_path, 'output')
        os.makedirs(output_path, exist_ok=True)

        # Get all file inputs and weights
        file_input = request.files['result_files[]']
        weight = request.form['weights[]']

        if not file_input.filename or not weight:
            raise ValueError("Both file and weight must be provided.")

        # Get the next available additional folder number
        existing_folders = [d for d in os.listdir(output_path) 
                          if os.path.isdir(os.path.join(output_path, d)) 
                          and d.startswith('additional_')]
        next_number = 1 if not existing_folders else max([int(d.split('_')[1]) for d in existing_folders]) + 1

        # Create new additional folder
        additional_folder_name = f'additional_{next_number}'
        additional_folder_path = os.path.join(output_path, additional_folder_name)
        os.makedirs(additional_folder_path, exist_ok=True)

        # Handle multiple files case
        if ',' in file_input.filename:
            # Split filenames and save all files
            filenames = file_input.filename.split(',')
            file_content = file_input.read()  # Read the file content once
            
            # Save all files
            for filename in filenames:
                filename = secure_filename(filename.strip())
                file_save_path = os.path.join(additional_folder_path, filename)
                with open(file_save_path, 'wb') as f:
                    f.write(file_content)  # Write the same content for each file
            
            # For input_list.txt, use only the first file after sorting
            sorted_filenames = sorted(filenames)
            first_filename = secure_filename(sorted_filenames[0].strip())
            first_file_path = os.path.join(additional_folder_path, first_filename)
            
            # Append to input_list.txt
            with open(os.path.join(output_path, 'input_list.txt'), 'a') as f:
                f.write(f"{first_file_path} {weight}\n")
        else:
            # Single file case
            filename = secure_filename(file_input.filename)
            file_save_path = os.path.join(additional_folder_path, filename)
            file_input.save(file_save_path)
            
            # Append to input_list.txt
            with open(os.path.join(output_path, 'input_list.txt'), 'a') as f:
                f.write(f"{file_save_path} {weight}\n")

        return jsonify(success=True)

    except Exception as e:
        print(f"Error in upload_files: {str(e)}")
        return jsonify(success=False, error=str(e))


@app.route('/upload_files', methods=['POST'])
def upload_files():
    try:
        working_directory = request.form.get('working_directory')
        if not working_directory:
            raise ValueError("Working directory not specified.")
        
        working_dir_path = os.path.join(UPLOAD_FOLDER, working_directory)
        output_path = os.path.join(working_dir_path, 'output')
        os.makedirs(output_path, exist_ok=True)

        # Get existing folder numbers
        existing_numbers = []
        for d in os.listdir(output_path):
            if os.path.isdir(os.path.join(output_path, d)) and d.startswith('additional_'):
                try:
                    num = int(d.split('_')[1])
                    existing_numbers.append(num)
                except (IndexError, ValueError):
                    continue
        
        next_number = 1 if not existing_numbers else max(existing_numbers) + 1
        
        # Create new additional folder
        additional_folder_name = f'additional_{next_number}'
        additional_folder_path = os.path.join(output_path, additional_folder_name)
        os.makedirs(additional_folder_path, exist_ok=True)

        # Get the file input and weight
        file_input = request.files['result_files[]']
        weight = request.form['weights[]']

        if not file_input.filename or not weight:
            raise ValueError("Both file and weight must be provided.")

        # Handle multiple files case
        if ',' in file_input.filename:
            # Save all files to the same folder
            filenames = file_input.filename.split(',')
            for filename in filenames:
                filename = secure_filename(filename.strip())
                file_save_path = os.path.join(additional_folder_path, filename)
                file_input.save(file_save_path)
            
            # For input_list.txt, use only the first file after sorting
            sorted_filenames = sorted(filenames)
            first_filename = secure_filename(sorted_filenames[0].strip())
            first_file_path = os.path.join(additional_folder_path, first_filename)
            
            # Append to input_list.txt
            with open(os.path.join(output_path, 'input_list.txt'), 'a') as f:
                f.write(f"{first_file_path} {weight}\n")
        else:
            # Single file case
            filename = secure_filename(file_input.filename)
            file_save_path = os.path.join(additional_folder_path, filename)
            file_input.save(file_save_path)
            
            # Append to input_list.txt
            with open(os.path.join(output_path, 'input_list.txt'), 'a') as f:
                f.write(f"{file_save_path} {weight}\n")

        return jsonify(success=True)

    except Exception as e:
        print(f"Error in upload_files: {str(e)}")
        return jsonify(success=False, error=str(e))

@app.route('/upload_files', methods=['POST'])
def upload_files():
    try:
        working_directory = request.form.get('working_directory')
        if not working_directory:
            raise ValueError("Working directory not specified.")
        
        working_dir_path = os.path.join(UPLOAD_FOLDER, working_directory)
        output_path = os.path.join(working_dir_path, 'output')
        os.makedirs(output_path, exist_ok=True)

        # Get existing folder numbers
        existing_numbers = []
        for d in os.listdir(output_path):
            if os.path.isdir(os.path.join(output_path, d)) and d.startswith('additional_'):
                try:
                    num = int(d.split('_')[1])
                    existing_numbers.append(num)
                except (IndexError, ValueError):
                    continue
        
        next_number = 1 if not existing_numbers else max(existing_numbers) + 1

        # Process each Upload Results to Merge section
        file_inputs = request.files.getlist('result_files[]')
        weight_inputs = request.form.getlist('weights[]')

        for file_input, weight in zip(file_inputs, weight_inputs):
            if not file_input.filename or not weight:
                continue

            # Create new additional folder
            additional_folder_name = f'additional_{next_number}'
            additional_folder_path = os.path.join(output_path, additional_folder_name)
            os.makedirs(additional_folder_path, exist_ok=True)

            if ',' in file_input.filename:
                # Multiple files case - save ALL files to the SAME folder
                filenames = file_input.filename.split(',')
                for filename in filenames:
                    filename = secure_filename(filename.strip())
                    file_save_path = os.path.join(additional_folder_path, filename)
                    file_input.save(file_save_path)
                
                # For input_list.txt, use only the first file after sorting
                sorted_filenames = sorted(filenames)
                first_filename = secure_filename(sorted_filenames[0].strip())
                first_file_path = os.path.join(additional_folder_path, first_filename)
                
                # Append to input_list.txt
                with open(os.path.join(output_path, 'input_list.txt'), 'a') as f:
                    f.write(f"{first_file_path} {weight}\n")
            else:
                # Single file case
                filename = secure_filename(file_input.filename)
                file_save_path = os.path.join(additional_folder_path, filename)
                file_input.save(file_save_path)
                
                # Append to input_list.txt
                with open(os.path.join(output_path, 'input_list.txt'), 'a') as f:
                    f.write(f"{file_save_path} {weight}\n")

            next_number += 1

        return jsonify(success=True)

    except Exception as e:
        print(f"Error in upload_files: {str(e)}")
        return jsonify(success=False, error=str(e))



@app.route('/upload_files', methods=['POST'])
def upload_files():
    try:
        working_directory = request.form.get('working_directory')
        if not working_directory:
            raise ValueError("Working directory not specified.")
        
        working_dir_path = os.path.join(UPLOAD_FOLDER, working_directory)
        output_path = os.path.join(working_dir_path, 'output')
        os.makedirs(output_path, exist_ok=True)

        # Get all file inputs and weights
        file_inputs = request.files.getlist('result_files[]')
        weight_inputs = request.form.getlist('weights[]')

        # Process each Upload Results to Merge section
        for file_input, weight in zip(file_inputs, weight_inputs):
            if not file_input.filename or not weight:
                continue

            # Get the next available additional folder number
            existing_folders = [d for d in os.listdir(output_path) 
                             if os.path.isdir(os.path.join(output_path, d)) 
                             and d.startswith('additional_')]
            next_number = len(existing_folders) + 1
            additional_folder_name = f'additional_{next_number}'
            additional_folder_path = os.path.join(output_path, additional_folder_name)
            os.makedirs(additional_folder_path, exist_ok=True)

            # Handle multiple files case
            if ',' in file_input.filename:
                # Save all files to the same additional folder
                filenames = file_input.filename.split(',')
                for filename in filenames:
                    filename = secure_filename(filename.strip())
                    file_save_path = os.path.join(additional_folder_path, filename)
                    file_input.save(file_save_path)
                
                # For input_list.txt, use only the first file after sorting
                sorted_filenames = sorted(filenames)
                first_filename = secure_filename(sorted_filenames[0].strip())
                first_file_path = os.path.join(additional_folder_path, first_filename)
                
                # Append to input_list.txt
                with open(os.path.join(output_path, 'input_list.txt'), 'a') as f:
                    f.write(f"{first_file_path} {weight}\n")
            else:
                # Single file case
                filename = secure_filename(file_input.filename)
                file_save_path = os.path.join(additional_folder_path, filename)
                file_input.save(file_save_path)
                
                # Append to input_list.txt
                with open(os.path.join(output_path, 'input_list.txt'), 'a') as f:
                    f.write(f"{file_save_path} {weight}\n")

        return jsonify(success=True)

    except Exception as e:
        print(f"Error in upload_files: {str(e)}")
        return jsonify(success=False, error=str(e))



@app.route('/upload_files', methods=['POST'])
def upload_files():
    try:
        working_directory = request.form.get('working_directory')
        if not working_directory:
            raise ValueError("Working directory not specified.")
        
        working_dir_path = os.path.join(UPLOAD_FOLDER, working_directory)
        output_path = os.path.join(working_dir_path, 'output')
        os.makedirs(output_path, exist_ok=True)

        # Get all file inputs and weights
        file_inputs = request.files.getlist('result_files[]')
        weight_inputs = request.form.getlist('weights[]')

        # Process each file input and weight pair
        for file_input, weight in zip(file_inputs, weight_inputs):
            if not file_input.filename or not weight:
                continue

            # Get the next available additional folder number
            existing_folders = [d for d in os.listdir(output_path) 
                             if os.path.isdir(os.path.join(output_path, d)) 
                             and d.startswith('additional_')]
            next_number = len(existing_folders) + 1
            additional_folder_name = f'additional_{next_number}'
            additional_folder_path = os.path.join(output_path, additional_folder_name)
            os.makedirs(additional_folder_path, exist_ok=True)

            # Handle multiple files case
            if ',' in file_input.filename:
                # Split filenames and save all files
                filenames = file_input.filename.split(',')
                for filename in filenames:
                    filename = secure_filename(filename.strip())
                    file_save_path = os.path.join(additional_folder_path, filename)
                    # Need to save the file content here
                    with open(file_save_path, 'wb') as f:
                        file_input.seek(0)
                        f.write(file_input.read())
            else:
                # Single file case
                filename = secure_filename(file_input.filename)
                file_save_path = os.path.join(additional_folder_path, filename)
                file_input.save(file_save_path)

        return jsonify(success=True)

    except Exception as e:
        print(f"Error in upload_files: {str(e)}")
        return jsonify(success=False, error=str(e))



@app.route('/upload_files', methods=['POST'])
def upload_files():
    try:
        working_directory = request.form.get('working_directory')
        if not working_directory:
            raise ValueError("Working directory not specified.")
        
        working_dir_path = os.path.join(UPLOAD_FOLDER, working_directory)
        output_path = os.path.join(working_dir_path, 'output')
        os.makedirs(output_path, exist_ok=True)

        # Get all file inputs and weights
        file_inputs = request.files.getlist('result_files[]')
        weight_inputs = request.form.getlist('weights[]')

        # Process each file input and weight pair
        for file_input, weight in zip(file_inputs, weight_inputs):
            if not file_input.filename or not weight:
                continue

            # Get the next available additional folder number
            existing_folders = [d for d in os.listdir(output_path) 
                             if os.path.isdir(os.path.join(output_path, d)) 
                             and d.startswith('additional_')]
            next_number = len(existing_folders) + 1
            additional_folder_name = f'additional_{next_number}'
            additional_folder_path = os.path.join(output_path, additional_folder_name)
            os.makedirs(additional_folder_path, exist_ok=True)

            if ',' in file_input.filename:  # Multiple files case
                # Save all files to the additional folder
                filenames = file_input.filename.split(',')
                for filename in filenames:
                    filename = secure_filename(filename.strip())
                    file_save_path = os.path.join(additional_folder_path, filename)
                    file_input.save(file_save_path)
                
                # For input_list.txt, use only the first file after sorting
                sorted_filenames = sorted(filenames)
                first_filename = secure_filename(sorted_filenames[0].strip())
                input_list_path = os.path.join(additional_folder_path, first_filename)
                
                # Append to input_list.txt
                with open(os.path.join(output_path, 'input_list.txt'), 'a') as f:
                    f.write(f"{input_list_path} {weight}\n")
            
            else:  # Single file case
                # Save the single file
                filename = secure_filename(file_input.filename)
                file_save_path = os.path.join(additional_folder_path, filename)
                file_input.save(file_save_path)
                
                # Append to input_list.txt
                with open(os.path.join(output_path, 'input_list.txt'), 'a') as f:
                    f.write(f"{file_save_path} {weight}\n")

        return jsonify(success=True)

    except Exception as e:
        print(f"Error in upload_files: {str(e)}")
        return jsonify(success=False, error=str(e))



@app.route('/upload_files', methods=['POST'])
def upload_files():
    try:
        working_directory = request.form.get('working_directory')
        if not working_directory:
            raise ValueError("Working directory not specified.")
        
        working_dir_path = os.path.join(UPLOAD_FOLDER, working_directory)
        output_path = os.path.join(working_dir_path, 'output')
        os.makedirs(output_path, exist_ok=True)

        # Get all file inputs and weights
        file_inputs = request.files.getlist('result_files[]')
        weight_inputs = request.form.getlist('weights[]')

        # Process each file input and weight pair
        for file_input, weight in zip(file_inputs, weight_inputs):
            if not file_input.filename or not weight:
                continue

            # For multiple files selected in a single input
            if ',' in file_input.filename:
                # Split and sort filenames
                filenames = file_input.filename.split(',')
                sorted_filenames = sorted(filenames)
                # Take only the first file after sorting
                filename = secure_filename(sorted_filenames[0])
            else:
                # Single file case
                filename = secure_filename(file_input.filename)

            # Get existing folder numbers
            existing_numbers = []
            for d in os.listdir(output_path):
                if os.path.isdir(os.path.join(output_path, d)) and d.startswith('additional_'):
                    try:
                        num = int(d.split('_')[1])
                        existing_numbers.append(num)
                    except (IndexError, ValueError):
                        continue

            # Calculate next folder number
            next_number = 1
            if existing_numbers:
                next_number = max(existing_numbers) + 1

            # Create new additional folder
            additional_folder_name = f'additional_{next_number}'
            additional_folder_path = os.path.join(output_path, additional_folder_name)
            os.makedirs(additional_folder_path, exist_ok=True)

            # Save the file
            file_save_path = os.path.join(additional_folder_path, filename)
            file_input.save(file_save_path)

            # Write to input_list.txt (append mode)
            with open(os.path.join(output_path, 'input_list.txt'), 'a') as f:
                f.write(f"{file_save_path} {weight}\n")

        return jsonify(success=True)

    except Exception as e:
        print(f"Error in upload_files: {str(e)}")
        return jsonify(success=False, error=str(e))


@app.route('/upload_files', methods=['POST'])
def upload_files():
    try:
        working_directory = request.form.get('working_directory')
        if not working_directory:
            raise ValueError("Working directory not specified.")
        
        working_dir_path = os.path.join(UPLOAD_FOLDER, working_directory)
        output_path = os.path.join(working_dir_path, 'output')
        os.makedirs(output_path, exist_ok=True)

        # Get all file inputs and weights
        file_inputs = request.files.getlist('result_files[]')
        weight_inputs = request.form.getlist('weights[]')

        # Process each file input and weight pair
        for file_input, weight in zip(file_inputs, weight_inputs):
            if not file_input.filename or not weight:
                continue

            # Get existing folder numbers and calculate next number
            existing_numbers = []
            for d in os.listdir(output_path):
                if os.path.isdir(os.path.join(output_path, d)) and d.startswith('additional_'):
                    try:
                        num = int(d.split('_')[1])
                        existing_numbers.append(num)
                    except (IndexError, ValueError):
                        continue
            
            next_number = 1
            if existing_numbers:
                next_number = max(existing_numbers) + 1

            # Create new additional folder
            additional_folder_name = f'additional_{next_number}'
            additional_folder_path = os.path.join(output_path, additional_folder_name)
            os.makedirs(additional_folder_path, exist_ok=True)

            # Handle the file(s)
            if ',' in file_input.filename:  # Multiple files were selected
                filenames = file_input.filename.split(',')
                sorted_filenames = sorted(filenames)
                selected_filename = sorted_filenames[0]  # Take first file after sorting
                filename = secure_filename(selected_filename)
            else:
                filename = secure_filename(file_input.filename)

            # Save the file
            file_save_path = os.path.join(additional_folder_path, filename)
            file_input.save(file_save_path)

            # Append to input_list.txt
            with open(os.path.join(output_path, 'input_list.txt'), 'a') as f:
                f.write(f"{file_save_path} {weight}\n")

        return jsonify(success=True)

    except Exception as e:
        print(f"Error in upload_files: {str(e)}")
        return jsonify(success=False, error=str(e))


@app.route('/upload_files', methods=['POST'])
def upload_files():
    try:
        working_directory = request.form.get('working_directory')
        if not working_directory:
            raise ValueError("Working directory not specified.")
        
        working_dir_path = os.path.join(UPLOAD_FOLDER, working_directory)
        output_path = os.path.join(working_dir_path, 'output')
        os.makedirs(output_path, exist_ok=True)

        # Get all file inputs and weights
        file_input = request.files['result_files[]']  # 단일 파일 입력 처리
        weight = request.form['weights[]']  # 단일 weight 처리

        if not file_input.filename or not weight:
            raise ValueError("Both file and weight must be provided.")

        # Get the next available additional folder number
        existing_folders = [d for d in os.listdir(output_path) 
                          if os.path.isdir(os.path.join(output_path, d)) 
                          and d.startswith('additional_')]
        next_number = len(existing_folders) + 1

        # Create new additional folder
        additional_folder_name = f'additional_{next_number}'
        additional_folder_path = os.path.join(output_path, additional_folder_name)
        os.makedirs(additional_folder_path, exist_ok=True)

        # Handle multiple files case
        if ',' in file_input.filename:
            # Split and sort filenames
            filenames = file_input.filename.split(',')
            sorted_filenames = sorted(filenames)
            selected_filename = sorted_filenames[0]  # Take first file after sorting
            filename = secure_filename(selected_filename)
        else:
            # Single file case
            filename = secure_filename(file_input.filename)

        # Save the file
        file_save_path = os.path.join(additional_folder_path, filename)
        file_input.save(file_save_path)

        # Write to input_list.txt
        with open(os.path.join(output_path, 'input_list.txt'), 'a') as f:
            f.write(f"{file_save_path} {weight}\n")

        return jsonify(success=True)

    except Exception as e:
        print(f"Error in upload_files: {str(e)}")  # 디버깅용 로그
        return jsonify(success=False, error=str(e))

@app.route('/upload_files', methods=['POST'])
def upload_files():
    try:
        working_directory = request.form.get('working_directory')
        if not working_directory:
            raise ValueError("Working directory not specified.")
        
        working_dir_path = os.path.join(UPLOAD_FOLDER, working_directory)
        output_path = os.path.join(working_dir_path, 'output')
        os.makedirs(output_path, exist_ok=True)

        # Get all file inputs and weights
        file_inputs = request.files.getlist('result_files[]')
        weight_inputs = request.form.getlist('weights[]')

        # Process each file input and weight pair independently
        for i, (file_input, weight) in enumerate(zip(file_inputs, weight_inputs)):
            if not file_input.filename or not weight:
                continue

            # Get the next available additional folder number
            existing_folders = [d for d in os.listdir(output_path) 
                             if os.path.isdir(os.path.join(output_path, d)) 
                             and d.startswith('additional_')]
            next_number = len(existing_folders) + 1

            # Create new additional folder
            additional_folder_name = f'additional_{next_number}'
            additional_folder_path = os.path.join(output_path, additional_folder_name)
            os.makedirs(additional_folder_path, exist_ok=True)

            # For multiple files in a single input
            if ',' in file_input.filename:
                filenames = file_input.filename.split(',')
                sorted_filenames = sorted(filenames)
                selected_filename = sorted_filenames[0]  # Take first file after sorting
                filename = secure_filename(selected_filename)
            else:
                filename = secure_filename(file_input.filename)

            # Save the file
            file_save_path = os.path.join(additional_folder_path, filename)
            file_input.save(file_save_path)

            # Write to input_list.txt
            with open(os.path.join(output_path, 'input_list.txt'), 'a') as f:
                f.write(f"{file_save_path} {weight}\n")

        return jsonify(success=True)

    except Exception as e:
        return jsonify(success=False, error=str(e))


@app.route('/upload_files', methods=['POST'])
def upload_files():
    try:
        working_directory = request.form.get('working_directory')
        if not working_directory:
            raise ValueError("Working directory not specified.")
        
        working_dir_path = os.path.join(UPLOAD_FOLDER, working_directory)
        output_path = os.path.join(working_dir_path, 'output')
        os.makedirs(output_path, exist_ok=True)

        # Get all file inputs and weights
        file_inputs = request.files.getlist('result_files[]')
        weight_inputs = request.form.getlist('weights[]')

        # Process each file input and weight pair
        for file_input, weight in zip(file_inputs, weight_inputs):
            if not file_input.filename or not weight:
                continue

            # Get current number of additional folders
            existing_folders = [d for d in os.listdir(output_path) 
                             if os.path.isdir(os.path.join(output_path, d)) 
                             and d.startswith('additional_')]
            next_number = len(existing_folders) + 1

            # Create new additional folder
            additional_folder_name = f'additional_{next_number}'
            additional_folder_path = os.path.join(output_path, additional_folder_name)
            os.makedirs(additional_folder_path, exist_ok=True)

            # Handle multiple files case
            if ',' in file_input.filename:
                # Split filenames and sort them
                filenames = file_input.filename.split(',')
                sorted_filenames = sorted(filenames)
                # Take only the first file after sorting
                filename = secure_filename(sorted_filenames[0])
                file_save_path = os.path.join(additional_folder_path, filename)
                # Save the first file
                file_input.save(file_save_path)
            else:
                # Single file case
                filename = secure_filename(file_input.filename)
                file_save_path = os.path.join(additional_folder_path, filename)
                file_input.save(file_save_path)

            # Append to input_list.txt
            with open(os.path.join(output_path, 'input_list.txt'), 'a') as f:
                f.write(f"{file_save_path} {weight}\n")

        return jsonify(success=True)

    except Exception as e:
        return jsonify(success=False, error=str(e))





@app.route('/upload_files', methods=['POST'])
def upload_files():
    try:
        working_directory = request.form.get('working_directory')
        if not working_directory:
            raise ValueError("Working directory not specified.")
        
        working_dir_path = os.path.join(UPLOAD_FOLDER, working_directory)
        output_path = os.path.join(working_dir_path, 'output')
        os.makedirs(output_path, exist_ok=True)

        # Get all file inputs and weights
        file_inputs = request.files.getlist('result_files[]')
        weight_inputs = request.form.getlist('weights[]')

        # Process each file input and weight pair
        for file_input, weight in zip(file_inputs, weight_inputs):
            if not file_input.filename or not weight:
                continue

            # Get list of existing additional folders
            existing_folders = [d for d in os.listdir(output_path) 
                             if os.path.isdir(os.path.join(output_path, d)) 
                             and d.startswith('additional_')]
            next_number = len(existing_folders) + 1

            # If multiple files selected, sort them
            if ',' in file_input.filename:
                filenames = file_input.filename.split(',')
                sorted_filenames = sorted(filenames)
                selected_filename = sorted_filenames[0]
            else:
                selected_filename = file_input.filename

            # Create new additional folder with unique number
            additional_folder_name = f'additional_{next_number}'
            additional_folder_path = os.path.join(output_path, additional_folder_name)
            os.makedirs(additional_folder_path, exist_ok=True)

            # Save the file
            filename = secure_filename(selected_filename)
            file_save_path = os.path.join(additional_folder_path, filename)
            
            # For multiple files case, need to find and save the correct file
            if ',' in file_input.filename:
                # Save the content of the first file
                file_input.save(file_save_path)
            else:
                file_input.save(file_save_path)

            # Append to input_list.txt
            with open(os.path.join(output_path, 'input_list.txt'), 'a') as f:
                f.write(f"{file_save_path} {weight}\n")

        return jsonify(success=True)

    except Exception as e:
        return jsonify(success=False, error=str(e))


@app.route('/upload_files', methods=['POST'])
def upload_files():
    try:
        working_directory = request.form.get('working_directory')
        if not working_directory:
            raise ValueError("Working directory not specified.")
        
        working_dir_path = os.path.join(UPLOAD_FOLDER, working_directory)
        output_path = os.path.join(working_dir_path, 'output')
        os.makedirs(output_path, exist_ok=True)

        # Get all file inputs and weights as separate groups
        file_inputs = request.files.getlist('result_files[]')
        weight_inputs = request.form.getlist('weights[]')

        # Process each file input and weight pair separately
        for file_input, weight in zip(file_inputs, weight_inputs):
            if not file_input.filename or not weight:
                continue

            # Handle multiple files in a single input
            if hasattr(file_input, 'getlist'):
                # Sort only within this group of files
                files = sorted(file_input.getlist(), key=lambda x: x.filename)
                if files:
                    file_to_process = files[0]  # Take the first file after sorting
            else:
                file_to_process = file_input

            # Create new additional folder
            existing_folders = [d for d in os.listdir(output_path) 
                             if os.path.isdir(os.path.join(output_path, d)) 
                             and d.startswith('additional_')]
            next_number = len(existing_folders) + 1
            additional_folder_name = f'additional_{next_number}'
            additional_folder_path = os.path.join(output_path, additional_folder_name)
            os.makedirs(additional_folder_path, exist_ok=True)

            # Save the file
            filename = secure_filename(file_to_process.filename)
            file_save_path = os.path.join(additional_folder_path, filename)
            file_to_process.save(file_save_path)

            # Append to input_list.txt
            with open(os.path.join(output_path, 'input_list.txt'), 'a') as f:
                f.write(f"{file_save_path} {weight}\n")

        return jsonify(success=True)

    except Exception as e:
        return jsonify(success=False, error=str(e))

@app.route('/upload_files', methods=['POST'])
def upload_files():
    try:
        working_directory = request.form.get('working_directory')
        if not working_directory:
            raise ValueError("Working directory not specified.")
        
        working_dir_path = os.path.join(UPLOAD_FOLDER, working_directory)
        output_path = os.path.join(working_dir_path, 'output')
        os.makedirs(output_path, exist_ok=True)

        # Get all file inputs and weights
        all_files = request.files.getlist('result_files[]')
        all_weights = request.form.getlist('weights[]')

        # Filter out empty file inputs
        file_weight_pairs = []
        for files, weight in zip(all_files, all_weights):
            if files.filename and weight:  # Only process if both file and weight exist
                if ',' in files.filename:  # Multiple files selected
                    filenames = files.filename.split(',')
                    sorted_filenames = sorted(filenames)
                    file_weight_pairs.append((files, weight))
                else:  # Single file selected
                    file_weight_pairs.append((files, weight))

        if not file_weight_pairs:
            raise ValueError("No valid file and weight pairs found.")

        # Process each valid file-weight pair
        for idx, (file, weight) in enumerate(file_weight_pairs):
            # Create new additional folder for each pair
            additional_folder_name = f'additional_{len(os.listdir(output_path)) + 1}'
            additional_folder_path = os.path.join(output_path, additional_folder_name)
            os.makedirs(additional_folder_path, exist_ok=True)

            if ',' in file.filename:  # Multiple files selected
                filenames = file.filename.split(',')
                sorted_filenames = sorted(filenames)
                filename = secure_filename(sorted_filenames[0])
            else:
                filename = secure_filename(file.filename)

            # Save the file
            file_save_path = os.path.join(additional_folder_path, filename)
            file.save(file_save_path)

            # Append to input_list.txt
            with open(os.path.join(output_path, 'input_list.txt'), 'a') as f:
                f.write(f"{file_save_path} {weight}\n")

        return jsonify(success=True)

    except Exception as e:
        return jsonify(success=False, error=str(e))
@app.route('/upload_files', methods=['POST'])
def upload_files():
    try:
        working_directory = request.form.get('working_directory')
        if not working_directory:
            raise ValueError("Working directory not specified.")
        
        working_dir_path = os.path.join(UPLOAD_FOLDER, working_directory)
        output_path = os.path.join(working_dir_path, 'output')
        os.makedirs(output_path, exist_ok=True)

        # Handle uploaded files for merging
        uploaded_files = request.files.getlist('result_files[]')
        weight_input = request.form.get('weights[]')

        if not uploaded_files:
            raise ValueError("No files uploaded.")
        
        if not weight_input:
            raise ValueError("No weight provided.")

        # If multiple files are uploaded, sort them and get the first one
        if len(uploaded_files) > 1:
            # Sort files by name
            sorted_files = sorted(uploaded_files, key=lambda x: secure_filename(x.filename))
            file_to_save = sorted_files[0]  # Take the first file after sorting
        else:
            file_to_save = uploaded_files[0]  # If only one file, use it directly

        # Create new additional folder
        additional_folder_name = f'additional_{len(os.listdir(output_path)) + 1}'
        additional_folder_path = os.path.join(output_path, additional_folder_name)
        os.makedirs(additional_folder_path, exist_ok=True)
        
        # Save the selected file
        filename = secure_filename(file_to_save.filename)
        file_save_path = os.path.join(additional_folder_path, filename)
        file_to_save.save(file_save_path)

        # Write to input_list.txt
        with open(os.path.join(output_path, 'input_list.txt'), 'a') as f:
            f.write(f"{file_save_path} {weight_input}\n")

        return jsonify(success=True)

    except Exception as e:
        return jsonify(success=False, error=str(e))
"""

@app.route('/save_input_list', methods=['POST'])
def save_input_list():
    try:
        data = request.get_json()
        working_directory = data.get('working_directory')
        input_list_content = data.get('input_list_content')
        
        if not working_directory:
            return jsonify(success=False, error="Working directory not specified."), 400
        
        input_list_path = os.path.join(UPLOAD_FOLDER, working_directory, 'output', 'input_list.txt')
        os.makedirs(os.path.dirname(input_list_path), exist_ok=True)

        with open(input_list_path, 'w') as f:
            f.write(input_list_content)

        return jsonify(success=True)

    except Exception as e:
        return jsonify(success=False, error=str(e)), 500

@app.route('/load_input_list', methods=['GET'])
def load_input_list():
    try:
        working_directory = request.args.get('working_directory')
        
        if not working_directory:
            return jsonify(success=False, error="Working directory not specified."), 400
        
        input_list_path = os.path.join(UPLOAD_FOLDER, working_directory, 'output', 'input_list.txt')

        if not os.path.exists(input_list_path):
            return jsonify(success=False, error="input_list.txt does not exist."), 404

        with open(input_list_path, 'r') as f:
            input_list_content = f.read()

        return jsonify(success=True, content=input_list_content)

    except Exception as e:
        return jsonify(success=False, error=str(e)), 500

@app.route('/download_merged_results/<working_directory>', methods=['GET'])
def download_merged_results(working_directory):
    try:
        merge_result_path = os.path.join(UPLOAD_FOLDER, working_directory, 'output', 'Merge', 'result.all.alleles')
        
        if not os.path.exists(merge_result_path):
            return jsonify(success=False, error="Merged result file does not exist."), 404
        
        return send_file(merge_result_path, as_attachment=True)

    except Exception as e:
        print(f"Error in /download_merged_results: {str(e)}")
        return jsonify(success=False, error="An error occurred while processing your request."), 500

# Run the Flask app
if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)

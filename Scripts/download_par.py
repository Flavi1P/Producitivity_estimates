import os
import requests
#python obdaac_download.py -v --filelist "C:/Users/flapet/OneDrive - NOC/Documents/IDAPro/lib/db_building/data/satellite/http_manifest.txt" --appkey 1526fffdc6e3c11287169a4ecf9bd3d22279fddd --uncompress --force --odir "C:/Users/flapet/OneDrive - NOC/Documents/IDAPro/lib/db_building/data/satellite/par2/"

# Directory to save the downloaded files
output_dir = "C:/Users/petit/Documents/Producitivity_estimates/Data/Raw/Satellite/PAR/"  # Replace with your desired directory
os.makedirs(output_dir, exist_ok=True)  # Create the directory if it doesn't exist

# List of URLs to download
urls = [
    "https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/requested_files_1.tar?h=ocdist408&p=/data1/60200244af799f6f/requested_files",
    # Add more URLs here
]

# Function to download a single file
def download_file(url, output_dir):
    filename = url.split("/")[-1]  # Extract the filename from the URL
    output_path = os.path.join(output_dir, filename)
    
    try:
        print(f"Downloading {filename}...")
        response = requests.get(url, stream=True)  # Stream the file
        response.raise_for_status()  # Raise an exception for HTTP errors
        
        with open(output_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):  # Download in chunks
                f.write(chunk)
        print(f"Downloaded: {output_path}")
    except requests.exceptions.RequestException as e:
        print(f"Failed to download {url}. Error: {e}")

# Download all files
for url in urls:
    download_file(url, output_dir)

"""
This script compresses the content of tmQM-RDF into .zip archives.

All the RDFS headers (data/tmQM-RDFS-*.ttl) are compressed into a single archive: tmQM-RDFS.zip.
All the RDF graphs (data/graphs/*.ttl) are compressed according to the first letter of their name into files of the form graphs/X.zip.
"""

# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))
    
# %% Header
import zipfile
import re

DATA = os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "data", "v2025dev")
RDFS_HEADER_PATTERN = "tmQM-RDFS-.\\.ttl"

# %% Utility functions

def compress_rdfs_headers(data_dir, rdfs_header_pattern):
    """
    This function compresses all the tmQM-RDFS headers into a single archive called tmQM-RDFS.zip
    (stored in [data_dir]/archive)
    
    - Arguments:
        - data_dir: the directory where the data is to be found
        - rdfs_header_pattern: a regexp specifying how to recognise RDFS headers
    """
    
    archive_dir = os.path.join(data_dir, "archive")
    if not os.path.exists(archive_dir):
        os.makedirs(archive_dir)
    
    # Locate RDFS header
    header_files = [x for x in os.listdir(data_dir) if re.match(rdfs_header_pattern, x) is not None]

    # Create zip archive
    with zipfile.ZipFile(os.path.join(archive_dir, "tmQM-RDFS.zip"), "w", compression = zipfile.ZIP_DEFLATED) as zipf:
        for hf in header_files:
            print(f"Compressing {hf}...")
            zipf.write(os.path.join(data_dir, hf), arcname = hf)
            
def compress_rdf_graphs(data_dir):
    """
    Compresses all the RDF graphs in tmQM-RDF into archives according to their initial letter.
    
    - Arguments:
        - data_dir: the directory where the data is to be found
    """
    
    archive_dir = os.path.join(data_dir, "archive", "graphs")
    if not os.path.exists(archive_dir):
        os.makedirs(archive_dir)
        
    # Organise graphs by initial
    graphs_dir = os.path.join(data_dir, "graphs")
    by_letter = {}
    
    for g in os.listdir(graphs_dir):
        if not g.endswith(".ttl"):
            continue
        
        if not g[0] in by_letter:
            by_letter[g[0]] = []
            
        by_letter[g[0]] += [g]
        
    # Create compressed archives
    for letter in by_letter:
        print(f"Compressing letter {letter}...")
        
        with zipfile.ZipFile(os.path.join(archive_dir, f"{letter}.zip"), "w", compression = zipfile.ZIP_DEFLATED) as zipf:
            for g in by_letter[letter]:
                zipf.write(os.path.join(graphs_dir, g), arcname = g)

# %% Main statement
if __name__ == "__main__":
    compress_rdfs_headers(DATA, RDFS_HEADER_PATTERN)
    compress_rdf_graphs(DATA)
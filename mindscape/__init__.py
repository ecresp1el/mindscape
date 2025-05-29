# MindScape package initialization

#imports the os module to interact with enviroment variables 
import os

# Enable debugging if the DEBUG environment variable is set
DEBUG = True and "DEBUG" in os.environ and os.environ["DEBUG"]

# Import version information
from mindscape.version import __version__, VERSION

# Print a message when the package is loaded
print(f"Loading MindScape {VERSION}...")
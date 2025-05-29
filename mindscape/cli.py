import click
from mindscape.create_project.new import create_new_project

@click.group()
def main():
    """MindScape Command-Line Interface"""
    pass

@main.command()
@click.argument("project")
@click.argument("experimenter")
@click.option("--working-directory", default=None, help="Base directory for the project.")
def create_project(project, experimenter, working_directory):
    """Create a new MindScape project."""
    create_new_project(project, experimenter, working_directory)
    print(f"✅ Project '{project}' created successfully!")

if __name__ == "__main__":
    main()
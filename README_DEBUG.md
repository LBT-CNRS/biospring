(README still under construction - 27/02/2025)

# Develop & Debug BioSpring in a devcontainer using VSCode

## Requirements:

- Act: https://github.com/nektos/act
- Docker Engine: https://docs.docker.com/engine/install/
- VSCode extensions:
  - GitHub Local Actions: https://marketplace.visualstudio.com/items?itemName=SanjulaGanepola.github-local-actions
  - Dev Containers: https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers

## Description of new files/folders in this branch

- .artifacts/ : Artifacts of the workflow that contains locally the FreeSASA and MDDriver sources in order to BioSpring to be built
- .devcontainer/
  - devcontainer.json : configuration file for the devcontainer using the biospring-debug image
- .github/workflows/local-debug.yml : workflow file to build the image of BioSpring in debug mode
- .actrc : configuration file for Act / Docker engine
- README_DEBUG.md : this file

## Build Docker Image for Local Debug in dev container

- Go to GitHub Local Actions > Workflows and run `Build Docker Image for Local Debug`
- Command Palette (CTRL+MAJ+P) > `Dev Containers: Reopen in Container` and wait until the C++ extensions finished to be installed
- Setup Debug config, see example below:

Example of `launch.json`

    {
        "version": "0.2.0",
        "configurations": [
            {
                "name": "(gdb) Launch",
                "type": "cppdbg",
                "request": "launch",
                "program": "/biospring/bin/biospring",
                "args": [
                    "-s", "/workspaces/biospring/example/2particles/model.nc",
                    "-c", "/workspaces/biospring/example/2particles/param.msp",
                    "--wait"
                ],
                "stopAtEntry": false,
                "cwd": "${fileDirname}",
                "environment": [],
                "externalConsole": false,
                "MIMode": "gdb",
                "setupCommands": [
                    {
                        "description": "Enable pretty-printing for gdb",
                        "text": "-enable-pretty-printing",
                        "ignoreFailures": true
                    },
                    {
                        "description": "Set Disassembly Flavor to Intel",
                        "text": "-gdb-set disassembly-flavor intel",
                        "ignoreFailures": true
                    }
                ]
            }
        ]
    }

- Run debug


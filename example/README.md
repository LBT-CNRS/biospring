# Usage

## Client and server on the same machine

  - ```bash
	make prep
	make run
	```
  - Check the line in the shell and note the `HOSTNAME` and `PORT`:
      ```bash
      MDDriver > Interactive MD bind to /[HOSTNAME]/[PORT]
      ```
  
  - In VMD: open `model.pdb`
  - In UnityMol:  open `model.nc` or `model.pdb`
  - Connect to the server using the `HOSTNAME` and `PORT`
  - CTRL + C to stop the server running the simulation
  - ```bash
	make clean
	```

## Client and server on different machines

- Step 1: Server side ⚙️

		make prep
		make expose_data


- Step 2: Client side 🖥️
  	
		make get_data
  
- Step 3: Server side ⚙️
  - CTRL + C to stop the server exposing the data
  - ```bash 
	make run
	```
  - Check the line in the shell and note the `HOSTNAME` and `PORT`:
      ```bash
      MDDriver > Interactive MD bind to /[HOSTNAME]/[PORT]
      ```
  
- Step 4: Client side 🖥️
  - In VMD: open `model.pdb`
  - In UnityMol:  open `model.nc` or `model.pdb`
  - Connect to the server using the `HOSTNAME` and `PORT`
  
- Step 5: Server side ⚙️
  - CTRL + C to stop the server running the simulation
  - ```bash
	make clean

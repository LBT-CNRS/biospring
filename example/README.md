# Usage:

## Step 1: Server side
- `make prep`
- `make expose_data`
  
## Step 2: Client side
- `make get_data`
  
## Step 3: Server side
- CTRL + C to stop the server exposing the data
- `make run` / `make run_now`
    - Check the line in the shell and note the `HOSTNAME` and `PORT`:
      ```bash
      MDDriver > Interactive MD bind to /[HOSTNAME]/[PORT]
      ```
  
## Step 4: Client side
- Open `model.nc` or `model.pdb` in VMD or UnityMol
- Connect to the server using the `HOSTNAME` and `PORT`
  
## Step 5: Server side
- CTRL + C to stop the server running the simulation
- `make clean`

import gurobipy as gp

def retrieve_wls_gurobi_license() -> dict:
    
    gurobi_license = {}
    
    # Open the gurobi.lic file in read mode
    with open('../gurobi.lic', 'r') as file:
        # Iterate over each line in the file
        for line in file:
            # Remove leading/trailing whitespace and newline characters
            line = line.split("=")
            
            # Use match and case to extract the desired values
            match line[0]:
                case "WLSACCESSID":
                    # Store the WLSACCESSID value
                    gurobi_license["WLSACCESSID"] = str(line[1].strip().replace("\n", ""))
                    wls_access_id = line[1]
                case "WLSSECRET":
                    # Store the WLSSECRET value
                    gurobi_license["WLSSECRET"] = str(line[1].strip().replace("\n", ""))
                    wls_secret = line[1]
                case "LICENSEID":
                    # Store the LICENSEID value
                    gurobi_license["LICENSEID"] = int(line[1].strip().replace("\n", ""))
                    license_id = line[1]
    
    return gurobi_license

def set_gurobi_key(gurobi_license: tuple):
    env = gp.Env(empty=True)
    env.setParam('WLSACCESSID', gurobi_license["WLSACCESSID"])
    env.setParam('WLSSECRET', gurobi_license["WLSSECRET"])
    env.setParam('LICENSEID', gurobi_license["LICENSEID"])
    env.start()

if __name__ == "__main__":
    # Initialize variables
    wls_access_id = None
    wls_secret = None
    license_id = None

    # Open the gurobi.lic file in read mode
    with open('gurobi.lic', 'r') as file:
        # Iterate over each line in the file
        for line in file:
            # Remove leading/trailing whitespace and newline characters
            line = line.split("=")
            
            # Use match and case to extract the desired values
            match line[0]:
                case "WLSACCESSID":
                    # Store the WLSACCESSID value
                    wls_access_id = line[1]
                case "WLSSECRET":
                    # Store the WLSSECRET value
                    wls_secret = line[1]
                case "LICENSEID":
                    # Store the LICENSEID value
                    license_id = int(line[1])

    # Print the extracted values
    print("WLSACCESSID:", type(wls_access_id))
    print("WLSSECRET:", type(wls_secret))
    print("LICENSEID:", type(license_id))

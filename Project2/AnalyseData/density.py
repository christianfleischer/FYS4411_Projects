import numpy as np
import matplotlib.pyplot as plt

def readData(filename):

    infile = open("Data/%s" %filename, 'r')
			
    positions = [[],[],[]]

    for line in infile:
        data = line.split()
        positions[0].append(float(data[0]))
        positions[1].append(float(data[1]))
        positions[2].append(float(data[2]))
    
    infile.close()

    return np.asarray(positions)
    
if __name__ == "__main__":
    positions = readData("positionsN10e6Int.dat")
    x = positions[0]
    y = positions[1]
    z = positions[2]
    
    r = np.sqrt(x**2 + y**2 + z**2)
    rMean = sum(r)/len(r)
    
    plt.hist(r, bins=100, normed=1)
    plt.axvline(rMean, color="r", label="Mean r = %f" %rMean)
    plt.legend()
    plt.title("N=10, MC cycles=1e6, Alpha=0.5, D=3")
    plt.xlabel(r"$r/a_{ho}$")
    plt.ylabel("Number of Samples (Normalized)")
    plt.show()

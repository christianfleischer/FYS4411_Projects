import sys
import numpy as np
import matplotlib.pyplot as plt

def readData(filename):
    energies = np.loadtxt("Data/%s" %filename)
    return energies
    
def blocking(energies, nBlocks, blockSize):
    
    meansOfBlocks = np.zeros(nBlocks)
    
    for i in range(nBlocks):
        energiesOfBlock = energies[i*blockSize:(i+1)*blockSize]
        meansOfBlocks[i] = sum(energiesOfBlock)/blockSize
    
    mean = sum(meansOfBlocks)/nBlocks
    mean2 = sum(meansOfBlocks**2)/nBlocks
    variance = mean2 - mean**2
    
    return mean, variance
    
if __name__ == "__main__":
    N = 1
    energies = readData("energiesN%i.dat" %N)
    
    deltaBlockSize = 100
    minBlockSize = 10
    maxBlockSize = 10000
    numberOfSizes = (maxBlockSize-minBlockSize)/deltaBlockSize + 1
    largestBlockSize = minBlockSize + (numberOfSizes-1)*deltaBlockSize #9910
    
    #blockSizes = np.zeros(numberOfSizes)
    blockSizes = np.linspace(minBlockSize, largestBlockSize, numberOfSizes).astype(int)
    blockAmounts = len(energies)/blockSizes#np.zeros(numberOfSizes)
    means = np.zeros(numberOfSizes)
    variances = np.zeros(numberOfSizes)
    
    for i in range(numberOfSizes):
        #blockSize = minBlockSize + i*deltaBlockSize
        #blockAmount = len(energies)/blockSize
        #mean, variance = blocking(energies, blockAmount, blockSize)
        mean, variance = blocking(energies, blockAmounts[i], blockSizes[i])
        means[i] = mean
        variances[i] = variance
    
    standardDeviation = np.sqrt(abs(variances)/(blockAmounts-1.))
    plt.plot(blockSizes, standardDeviation)
    plt.xlabel("Block Size")
    plt.ylabel(r"Standard Deviation $\sigma$")
    plt.title("N=%i" %N)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.show()
    


import cwt.cwt as cwt
import cwt.correctbounds as correctbounds 
import pylab as pl



def main():
    #visualize = False
    visualize = True
    minimumIntensity = 100.0




    #f = open('testing_eic.txt','r')
    f = open('testing_eic2.txt','r')
    d = pl.genfromtxt(f)

    intensities = d[0,:]
    rt = d[1,:]


    if visualize:
        pl.plot(rt,intensities)
        pl.show()

    cwtObject = cwt.ContinuousWaveletTransform(1,10,1)
    cwtObject.setVisualize(visualize )
    cwtObject.setSignal(intensities)
    cwtObject.setX(rt)
    cwtObject.buildRidgelines()
    cwtObject.filterRidgelines()

    newPeaks = cwtObject.findBoundries();


    for i in range(len(newPeaks[0,:])):

        leftBound = int(newPeaks[0,i])
        rightBound = int(newPeaks[1,i])

        if visualize:
            pl.plot(rt,intensities)
            pl.fill_between(rt[leftBound:rightBound+1],intensities[leftBound:rightBound+1],facecolor='r')
            pl.title('peak before boundry correction')
            pl.show()

        leftBound = correctbounds.fixLeftBoundry(intensities,leftBound )
        rightBound = correctbounds.fixRightBoundry(intensities,rightBound )

        leftBound,rightBound = correctbounds.cropZerosFromEdges(intensities,leftBound,rightBound)

        if visualize:
            pl.plot(rt,intensities)
            pl.fill_between(rt[leftBound:rightBound+1],intensities[leftBound:rightBound+1],facecolor='r')
            pl.title('peak after boundry correction')
            pl.show()

        print "left bound = " + str(leftBound)
        print "right bound = " + str(rightBound)

        if max(intensities[leftBound:rightBound])<minimumIntensity:
            continue

if __name__ == "__main__":
    main()

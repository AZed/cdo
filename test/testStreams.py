from cdo import *
def createInputFiles(nfiles,resolution,cdo,cdoOptions):
  timeDefaults  = '-settunits,days -settaxis,20010101,12:00:00'
  randomEnlarge = lambda  resolution: '-mul -random,'+resolution+' -enlarge,'+resolution
  fileList      = {}
  cdo.debug     = False

  FileKeys = ['1d','2d','3d','1d_withTime','2d_withTime','3d_withTime']
  for i,key in enumerate(FileKeys):
    fileList[key] = []
    if '1d' == key:
      # vertical levels only
      for j in xrange(0,nfiles):
        fileList[key].append(cdo.stdatm(0,10,20,50, options = cdoOptions))

    elif '2d' == key:
      for j in xrange(0,nfiles):
	fileList[key].append(cdo.setcode(130,input = "-setname,'T' -random,"+resolution, options = cdoOptions))

    elif '3d' == key:
      for j in xrange(0,nfiles):
	fileList[key].append(cdo.enlarge(resolution, input = '-stdatm,0,10,20,50', options = cdoOptions))

    elif '1d_withTime' == key:
      # temporal axis only
      for j in xrange(0,nfiles):
	fileList[key].append(cdo.setcode(130,input = cdo.setname('T',input = timeDefaults+' -for,1,400', options = cdoOptions)))

    elif '2d_withTime' == key:
      for j in xrange(0,nfiles):
	fileList[key].append(cdo.setcode(130,input = cdo.setname('T',input = randomEnlarge(resolution) +' ' +timeDefaults+ ' -for,1,400', options = cdoOptions)))

    elif '3d_withTime' == key:
      for j in xrange(0,nfiles):
	fileList[key].append(cdo.enlarge(resolution, input = '-mul '+timeDefaults+ ' -stdatm,0,10,20 -for,1,400', options = cdoOptions))

  return fileList


if __name__ == '__main__':
  cdo = Cdo()
  a = createInputFiles(5,'r36x18',cdo,'-f nc')
  for i,value in enumerate(a):
    files = a[value]
    for f in files:
      print(' '.join([value,f,os.path.exists(f).__str__()]))

import FWCore.ParameterSet.Config as cms

def restrictInput( fileNames, numberOfFiles, firstFile=0) :

    tmpFiles =  fileNames
    files = cms.untracked.vstring()
    i = 0
    for file in tmpFiles:
        if i>=firstFile:
            files.append( file )
        if len(files)== numberOfFiles:
            break
        i+=1
    print 'restricting input to:'
    print files 
    return files

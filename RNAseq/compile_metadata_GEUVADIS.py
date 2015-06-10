from optparse import OptionParser
import os
import math
from time import gmtime, strftime

'''
Author:Oana Ursu
'''

def main():
    parser=OptionParser()
    
    parser.add_option('--out',dest='out',help='Out')
    geuvadis_meta='/srv/gs1/projects/snyder/jzaugg/histoneQTL/RNAseq/data/E-GEUV-1.sdrf.txt'
    our_study='/srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment'
    pickrell='/srv/gs1/projects/snyder/jzaugg/histoneQTL/RNAseq/data/list_lanes_pickrell_2010_nature'
    opts,args=parser.parse_args()
    
    our_people=set()
    for line in open(our_study,'r').readlines():
        our_people.add(line.strip().split('\t')[0])
    print our_people
                   
    geu1=set()
    for line in open(geuvadis_meta,'r').readlines():
        items=line.strip().split('\t')
        geu1.add(items[0])
    print geu1

    print 'in common'
    print geu1.intersection(our_people)
    print str(len(geu1.intersection(our_people)))

    residual_individuals=our_people.difference(geu1)
    print 'residual individuals'
    print residual_individuals

    pickrell_people=set()
    print 'How many residual individuals in Degner?'
    for line in open(pickrell,'r').readlines():
        pickrell_people.add(line.strip().split()[3])
    print 'pickrell'
    print pickrell_people
    print 'gets these residuals'
    print len(residual_individuals.intersection(pickrell_people))

    print 'pickrell gets these individuals from our study'
    print len(our_people.intersection(pickrell_people))
    print our_people.intersection(pickrell_people)

    print 'geuvadis gets these from our study'
    print our_people.intersection(geu1)

    print 'overlap geuvadis and pickrell'
    print len(geu1.intersection(pickrell_people))
    print geu1.intersection(pickrell_people)

    geu1_and_our=geu1.intersection(our_people)
    print 'geuvadis and us'
    print len(geu1_and_our)
    print (geu1_and_our)
    print 'in geuvadis and our study but not pickrell'
    print len(geu1_and_our.difference(pickrell_people))
    print geu1_and_our.difference(pickrell_people)



main()

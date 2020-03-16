# -*- coding: utf-8 -*-
"""
Created on Fri Jul  6 12:31:03 2018

@author: sosarica
"""
import sys

nbx=[]
labels=[]

filename=sys.argv[1]
print (filename)



filein=open(filename)
print ("Reading file"+str(filein))
n=0
current_fill=0

for line in filein.readlines():
     n+=1
     if n<=2:
         labels.append(line)
     else:
         items= line.split(",")            
         if current_fill!=int(items[0]):
             nbx.append((int(items[0]),int(items[7])))
             current_fill=int(items[0])

filename=filename.split(".")[0]    
fileout=open(filename+"_perFill.csv","w")
    
for dat in nbx:
    fileout.write(str(int(dat[0]))+","+str(dat[1])+ "\n")
        
        
print ("This file contains: ",labels)
print (nbx)

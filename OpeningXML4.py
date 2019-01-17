########################################
##      Most Recent Code to Help      ##                         
##      Drew Open His XML Files       ##                        
########################################



from xml.dom import minidom
import numpy as np
from tqdm import tqdm

filetags = ["0p5","1p0","1p5","2p0","2p5","3p0","3p5","4p0","4p5", "5p0"]


for ct in tqdm(range(len(filetags))):

    currenttag = ct
    
    filename = 'C:/Users/PSU_Telemetry/Desktop/XML Stuff/drew_polar_scan/drew_polar_scan_gamma' + filetags[currenttag] + '_results.xml'
    file = open(filename,'r')
    doc = minidom.parse(file)
    file.close()
    
    string = doc.toxml()
    
    file = open('C:/Users/PSU_Telemetry/Desktop/XML Stuff/xmloutput.txt','w')
    file.write(doc.toprettyxml())
    file.close()
    
    numSteps = string.count("Step_number")
    
    anglist = []
    orders = [0] * numSteps
    efficiencies = [0] * numSteps
    
    for i in range(numSteps):
    
        #jump to scanning parameter value
        string = string[string.find("Scanning_parameter_value")+26:]
        
        #isolate scanning parameter value
        paramValue = float(string[:string.find('"')])
        
        anglist.append(paramValue)
        
        orderlist = []
        efflist = []
        
        partstring = string[:string.find("</Far_field>")]
        numOrders = partstring.count("Azimuth_angle")
        
        for j in range(numOrders):
        
            #jump to efficiency
            string = string[string.find("Efficiency_TE=")+15:]
            #isolate efficiency
            efficiency = float(string[:string.find("Exist_For")-2])
            efflist.append(efficiency)
            
            #jump to order
            string = string[string.find("Order_number")+14:]
            order = float(string[:string.find("Phase_TE")-2])
            orderlist.append(order)
        
        orders[i] = orderlist
        efficiencies[i] = efflist
    
    outputfilename = 'C:/Users/PSU_Telemetry/Desktop/XML Stuff/XML Files 4/drew_polar_scan_gamma' + filetags[currenttag] + '_results TEXT.txt'
    file = open(outputfilename,"w")
    
    file.write("Angles,Orders,Efficiencies")
    file.write("\n")
        
    for i in range(len(anglist)):
        file.write(str(anglist[i]) + ",")
        for item in orders[i]:
            file.write(str(item) + ",")
        for item in efficiencies[i]:
            file.write(str(item) + ",")
        file.write("\n")
    file.close()
    
    # Write 5th,6th,7th orders only
    outputfilename = 'C:/Users/PSU_Telemetry/Desktop/XML Stuff/XML Files 4/drew_polar_scan_gamma' + filetags[currenttag] + '_results 4th_5th orders only.txt'
    file = open(outputfilename,"w")
    
    file.write("Angles,Orders,Efficiencies")
    file.write("\n")
        
    for i in range(len(anglist)):
        
        arr = np.array(orders[i])
        wanted = np.logical_and((arr > 3),(arr < 6))
        wantedtemp = np.logical_and((arr == -4),(arr == -5))
        wanted = np.logical_or(wanted,wantedtemp)
        
        if ((4 in orders[i]) or (5 in orders[i])):
            file.write(str(anglist[i]) + ",")
        for item in np.array(orders[i])[wanted]:
            file.write(str(item) + ",")
        for item in np.array(efficiencies[i])[wanted]:
            file.write(str(item) + ",")
        if ((4 in orders[i]) or (5 in orders[i])):
            file.write("\n")
        
    file.close()





















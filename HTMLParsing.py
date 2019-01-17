import os
import io
papers = os.listdir("C:/Users/PSU_Telemetry/Desktop/Website Papers")

output = ""

data_in = open('titles.txt', 'r')
mid = data_in.read()
data = mid.split()
data_in.close()

data = "".join(data)
pieces = data.split("Filename:")
titles = []
for p in pieces:
    titles.append(p.split("Edit")[0])
titles = titles[1:]

numbers = []
for i in range(len(papers)):
    numbers.append(papers[i][:papers[i].find("-")])
    
file = open('Publication Order.txt','r')
order = file.read()
order = order.split("\n")
file.close()

#<a href="https://sites.psu.edu/mcentaffergroup/files/2018/09/1-GratingDesignForTheWaterRecoveryXRayRocket-1ksbtr1.pdf">Title Here </a>

output = []

for i in range(len(order)-4):
    
    if (order[i].find("@") != -1):
        text = order[i][order[i].find("@")+1:]
        output.append("<h2>" + text + "</h2><br>")
        continue
        
    papertitle = order[i]
    papernumber = papertitle[:papertitle.find("-")]
    paperindex = numbers.index(papernumber)
    title = titles[paperindex]
    paper = papers[paperindex][papers[paperindex].find("-")+1:-4]
    
    output.append('<a href="https://sites.psu.edu/mcentaffergroup/files/2018/09/' + title + '">')
    output.append(paper + "</a><br>")
    
with io.open("HTMLFile.txt","w",encoding='utf-8') as f:
    
    for o in output:
        f.write(o)
        f.write("\n")
    
import numpy as np

#perm = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49]
perm = [2,5,6,7,8,9,16,17,10,12,13,40,41,42,19,24,25,27,28,29,30,31,33,0,1,3,4,11,14,15,18,20,21,22,23,26,32,34,35,36,37,38,44,45,46,47,48,39,43,49]
flag = [1,1,1,1,1,-1,1,1,1,-1,1,1,1,-1,1,1,1,1,1,1,1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1,1,1,1,1,1,-1]
file = open("Ecoli.fas","w")
total = ""
    
for i in range(50):
    file2 = open("Ecoli_"+str(perm[i]+1)+".txt","r")
    text = file2.read().replace('\n','')
    if flag[perm[i]] == -1:
        text = text.replace("A","o")
        text = text.replace("T","A")
        text = text.replace("o","T")
        text = text.replace("a","o")
        text = text.replace("t","a")
        text = text.replace("o","t")
        text = text.replace("C","o")
        text = text.replace("G","C")
        text = text.replace("o","G")
        text = text.replace("c","o")
        text = text.replace("g","c")
        text = text.replace("o","g")
        total += text[::-1]
    else:
        total += text

    norm = "N"*(6000000-len(text))
    print(6000000-len(text))
    total += norm
    print(len(total))
    print(perm[i], flag[perm[i]])
    file2.close()

print(len(total))
file.write(total)
file.close()
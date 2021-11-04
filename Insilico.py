import numpy as np
import random
import time
import os
import keyboard
import threading

THRESHOLD = 3000

def macro():
    while True:
        keyboard.press_and_release('\n')
        time.sleep(0.5)
        
def zipfian(s,N):
    temp0 = np.array(range(1,N+1))
    temp0 = np.sum(1/temp0**s)
    temp = random.random() * temp0

    for i in range(N):
        temp2 = 1 / ((i + 1) ** s)

        if temp < temp2:
            return i+1
        else:
            temp = temp - temp2

    return 0

def ReadSeq(fname):
    file = open(fname,'r')
    SeqTemp = list(file.readline().replace("\n",""))
    Seq = np.zeros(np.size(SeqTemp))
    for _ in range(np.size(SeqTemp)):
        Seq[_] = (SeqTemp[_] == 'A') + 2 * (SeqTemp[_] == 'C') + 3 * (SeqTemp[_] == 'G') + 4 * (SeqTemp[_] == 'T') - 1
    file.close()
    
    return Seq

def Setparam(name, i, j):
    file = open("param.txt","w")

    file.write("INPUT_FILE	"+name+"	// input file name (name size limit: 256)\n")
    file.write("WORD_SIZE	14		// WORD size (W, supported to 16 )\n")
    file.write("ALLOW_SIZE	1		// allowable size per word (m, suport 0 or 1)\n")
    file.write("SPACE_SIZE	2		// SPACE size (SP)\n")
    file.write("MIN_SEED_LEN	28		// minimum length to be a seed (L)\n")
    file.write("SCORE_MAT	1		// match score (integer)\n")
    file.write("SCORE_MIS	-2		// mismatch score (integer)\n")
    file.write("SCORE_THR	-10		// score threshold (in ungapped extension)\n")
    file.write("GREEDY_X	30		// value of X in greedy algorithm (in gapped extension)\n")
    file.write("GREEDY_MIN_L	-10240		// maximum value of insertion\n")
    file.write("GREEDY_MAX_U	10240		// maximum value of deletion\n")
    file.write("WD_SIZE		20		// window size for filtering\n")
    file.write("T_THR		0.6		// threshold for filtering\n")

    file.close()
    
    os.system("mkdir Ecoli\\Ecoli_"+str(i)+"_"+str(j))
    
    file = open("param2.txt","w")

    file.write("REM_FILE	result.rem	 	// rem file\n")
    file.write("MAX_DIFF	20000		// 최대 차이값\n")
    file.write("SEQ_LINE	60		//  seq를 보여주는 라인 길이\n")
    file.write("ALIGN_PATH	Ecoli\\Ecoli_"+str(i)+"_"+str(j)+"		// alignment 결과파일을 생성할 디렉토리\n")
    file.write("ALIGN_SIZE	6000000		// alignment 결과파일 정리를 위한 기본 지역 크기\n")
    file.write("ALIGN_LENGTH	1000000		// alignment 결과파일 정리를 위한 기본 alignment 길이\n")
    file.write("ALIGN_MAX	1000000		// alignment 결과파일 정리를 위한 최대 alignment 길이\n")
    
    file.close()


def Insilico(seq, snp, indel, maxI, num, name):
    l_seqs = len(seq)
    seqs = []
    for i in range(num):
        past = time.time()
        seqtemp = np.mod(seq + (np.random.rand(l_seqs)  < snp)*np.random.randint(4, size=l_seqs),4)
        indelonoff = (np.random.rand(l_seqs) < indel)
        indeltotal = np.sum(indelonoff)
        indelindex = np.where(indelonoff == 1)
        deleteuntil = 0
        #print('total indels', indeltotal)
        
        now = time.time()
        #print('generate random array', now-past)
        past = time.time()
        indeltemp = []
        for kk in range(indeltotal):
            #if kk%10000 == 0:
                #print(i, 0, kk)
            if indelindex[0][kk] > deleteuntil:
                indellen = zipfian(1.6,maxI)
                ranval = np.random.rand()
                if ranval < 1/2:
                    indeltemp.append(np.random.randint(4, size=indellen))
                else:
                    indeltemp.append(-indellen)
                    deleteuntil = indelindex[0][kk] + indellen
        
        #print(indeltemp)
        now = time.time()
        #print('indel generation', now-past)
        past = time.time()
        lastindex = 0
        file = open(name,'w')
        for kk in range(len(indeltemp)):
            #if kk%10000 == 0:
                #print(i, 1, kk)
            if np.sum(indeltemp[kk]) < 0:
                seqtemp3 = seqtemp[lastindex:indelindex[0][kk]]
                lastindex = indelindex[0][kk]-indeltemp[kk]
                seq2 = ''
                for j in range(len(seqtemp3)):
                    if seqtemp3[j] == 0:
                        seq2 = seq2 + 'A'
                    elif seqtemp3[j] == 1:
                        seq2 = seq2 + 'C'
                    elif seqtemp3[j] == 2:
                        seq2 = seq2 + 'G'
                    elif seqtemp3[j] == 3:
                        seq2 = seq2 + 'T'
                file.write(seq2)
            else:
                seqtemp3 = np.append(seqtemp[lastindex:indelindex[0][kk]],indeltemp[kk])
                lastindex = indelindex[0][kk]
                seq2 = ''
                for j in range(len(seqtemp3)):
                    if seqtemp3[j] == 0:
                        seq2 = seq2 + 'A'
                    elif seqtemp3[j] == 1:
                        seq2 = seq2 + 'C'
                    elif seqtemp3[j] == 2:
                        seq2 = seq2 + 'G'
                    elif seqtemp3[j] == 3:
                        seq2 = seq2 + 'T'
                file.write(seq2)
                
        seqtemp3 = seqtemp[lastindex:l_seqs]
        seq2 = ''
        for j in range(len(seqtemp3)):
            if seqtemp3[j] == 0:
                seq2 = seq2 + 'A'
            elif seqtemp3[j] == 1:
                seq2 = seq2 + 'C'
            elif seqtemp3[j] == 2:
                seq2 = seq2 + 'G'
            elif seqtemp3[j] == 3:
                seq2 = seq2 + 'T'
        file.write(seq2)
        
        file.close()       
        now = time.time()
        #print('write sequence', now-past)

def AlignProcess(i, j):
    fname = "D:\\Dropbox\\Development\\REminerII\\Processed2\\Ecoli_"+str(i)+"_"+str(j)+".csv"
    dir1 = "D:\\Dropbox\\Development\\REminerII\\Ecoli_"+str(i)+"_"+str(j)+"\\Forward"
    dir2 = "D:\\Dropbox\\Development\\REminerII\\Ecoli_"+str(i)+"_"+str(j)+"\\Reverse"
    
    fname1 = os.listdir(dir1)
    fname2 = os.listdir(dir2)
    startlist = []
    endlist = []
    scorelist = []
    ACGTlist = []
    
    for name in fname1:
        fnamet = os.listdir(dir1+"\\"+name)
        for name2 in fnamet:
            fnamett = os.listdir(dir1+"\\"+name+"\\"+name2)
            for name3 in fnamett:
                if "align" in name3:
                    info = name3.split("(")[1].split(")")[0]
                    start = info.split("_")[0]
                    end = info.split("_")[1]
                    score = int(info.split("_")[2].split(".")[2])
                    if score > THRESHOLD:
                        file = open(dir1+"\\"+name+"\\"+name2+"\\"+name3,'r')
                        data = file.read()
                        ACGT = [data.count("A"),data.count("C"),data.count("G"),data.count("T")]
                        file.close()
                        print(start, end, score)
                        startlist.append(start.split("."))
                        endlist.append(end.split("."))
                        scorelist.append(score)
                        ACGTlist.append(ACGT)
                        
        fnamet = os.listdir(dir2+"\\"+name)
        for name2 in fnamet:
            fnamett = os.listdir(dir2+"\\"+name+"\\"+name2)
            for name3 in fnamett:
                if "align" in name3:
                    info = name3.split("(")[1].split(")")[0]
                    start = info.split("_")[0]
                    end = info.split("_")[1]
                    score = int(info.split("_")[2].split(".")[2])
                    if score > THRESHOLD:
                        file = open(dir2+"\\"+name+"\\"+name2+"\\"+name3,'r')
                        data = file.read()
                        ACGT = [data.count("A"),data.count("C"),data.count("G"),data.count("T")]
                        file.close()
                        print(start, end, score)
                        startlist.append(start.split("."))
                        endlist.append(end.split("."))
                        scorelist.append(-score)
                        ACGTlist.append(ACGT)
                        
    alignstart = []
        
    for _ in range(len(startlist)):
        alignstart.append([int(startlist[_][0]) , int(startlist[_][1])])
        
    alignend = []
        
    for _ in range(len(startlist)):
        alignend.append([int(endlist[_][0]), int(endlist[_][1])])
            
    alignscore = []
        
    for _ in range(len(startlist)):
        alignscore.append(scorelist[_])
            
    alignACGT = []
        
    for _ in range(len(startlist)):            
        alignACGT.append(ACGTlist[_])
    
    file = open(fname,'w')
    for k in range(len(alignscore)):
        xdiff = alignend[k][0] - alignstart[k][0]
        ydiff = alignend[k][1] - alignstart[k][1]
        diff = np.abs(xdiff) + np.abs(ydiff)
        file.write(str(alignstart[k][0]/6000000)+","+str(alignstart[k][1]/6000000)+","+str(alignend[k][0]/6000000)+","+str(alignend[k][1]/6000000)+","+str(2*alignscore[k]/diff)+","+str(alignACGT[k][0]/diff)+","+str(alignACGT[k][1]/diff)+","+str(alignACGT[k][2]/diff)+","+str(alignACGT[k][3]/diff)+"\n")
        
    file.close()
    

Macro = threading.Thread(target=macro)
Macro.start()

for i in range(50):
    past = time.time()
    seq = ReadSeq('Ecoli_'+str(i+1)+'.txt')
    
    now = time.time()
    past = now
    print('Read sequence', now-past)
    for j in range(10):
        Insilico(seq,0.05,0.01,100,1,"D:\\Dropbox\\Development\\REminerII\\sequence\\Ecoli_"+str(i)+"_"+str(j)+".txt")
        Setparam("D:\\Dropbox\\Development\\REminerII\\sequence\\Ecoli_"+str(i)+"_"+str(j)+".txt", i, j)
        os.system("REMINER2_VS2008.exe")
        time.sleep(0.1)
        os.system("Release\\AlignmentGenerator.exe")
        time.sleep(0.1)
        AlignProcess(i,j)
        
    now = time.time()
    print('Process sequence', now-past)
    
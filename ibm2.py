#Zachary Hinz, Abby Smith, Riki Imai, Maya Singapuri
import random


class ibm2:


 def __init__(self, e_sent, f_sent, iterations, threshhold, alignmentsflag):
     engFile, forFile= open(e_sent, "r"), open(f_sent, "r")


     # P(f|e), alignment_dict
     # P(e, f)/count(e, f), count (e)
     pfe, alignment = {}, {}
     engSent, forSent, fileLength = [], [], 0
  
     # read english File into arraylist
     for line in engFile:
         engSent.append("NULL "+line)
         #engSent.append(line)
         fileLength += 1


     # read foreign file into arraylist
     for line in forFile:
         forSent.append(line)


     engFile.close()
     forFile.close()


     engVocab = [] #added to make IBM model 2 easier
     forVocab = []
     # set up pfe tables, assigning them a value of 0.01
     for i in range(fileLength):
         for e in list(set(engSent[i].split())):
             engVocab.append(e)
             for f in list(set(forSent[i].split())):
                 forVocab.append(f)   
                 if(e not in pfe):
                     pfe[e] = {f:0.01}
                 elif(f not in pfe[e]):
                     pfe[e][f] = 0.01
                
     print("PFE Table Created")


     # complete alignment alg for given iterations
     for times in range(iterations):
         # instantiate the count tables
         ce = {}
         cef = {}
         for i in range(fileLength):
             # iterate through each pairing
             if(i%1000==0):
                 print("Line Num: "+str(i))
             engWords = list(set(engSent[i].split()))
             for e in engWords:
                 # add word to ce and cef tables as needed
                 if e not in ce:
                     ce[e] = 0
                     cef[e] = {}
                 forWords = list(set(forSent[i].split()))
                 cache = {}
                 for f in forWords:
                     # add word to cef[e] as needed
                     if f not in cef[e]:
                         cef[e][f] = 0
                     # cache to reduce calculations
                     if f in cache:
                         denom = cache[f]
                     else:
                         denom = 0
                         for e_s in engWords:
                             denom+=pfe[e_s][f]
                         cache[f] = denom
                  
                     pFtoE = pfe[e][f]/denom
                     cef[e][f] += pFtoE
                     ce[e] += pFtoE 


         # update p(f->e) after setup        
         for eng in cef:
             for foreign in cef.get(eng):
                 pfe.get(eng)[foreign] = (cef.get(eng)[foreign])/(ce[eng])




         print("Iteration: "+str(times))


      # ***** IBM MODEL 2 *****


      #carry over p(e|f) from model 1 (TODO is this the same as p(f|e??)

    #TODO a while loop for number of iterations?? if so how many

     print(pfe)
     print("start IMB 2")

     # line 2 of pseudocode -- initalize a(i|j, le, lf) to 1/(lf+1) for all i, j, le, lf
     alignModel2 = {} # goes [f][e][le][lf]
     for i in range(0, len(engSent)):
         lf = len(forSent[i].split())
         le = len(engSent[i].split())
         value = 1 / (lf)
         for f in range(0, lf):
             for e in range(0, le):
                 j = forSent[i].split()[f]
                 k = engSent[i].split()[e]

                 if(j not in alignModel2):
                     alignModel2[j] = {}
                 if(k not in alignModel2[j]):
                     alignModel2[j][k] = {}
                 if(le not in alignModel2[j][k]):
                     alignModel2[j][k][le] = {}
                 if(lf not in alignModel2[j][k][le]):
                     alignModel2[j][k][le][lf] = {}
                 alignModel2[j][k][le][lf] = value


     for times in range(iterations):

        countAlign = {}  #4 things, line 7
        totalAlign = {} #3 things, line 8
        for i in range(0, len(engSent)):
            lf = len(forSent[i].split())
            le = len(engSent[i].split())
            for f in range(0, lf):
                for e in range(0, le):
                    j = forSent[i].split()[f]
                    k = engSent[i].split()[e]
                    if(j not in countAlign):
                        countAlign[j] = {}
                    if(k not in countAlign[j]):
                        countAlign[j][k] = {}
                    if(le not in countAlign[j][k]):
                        countAlign[j][k][le] = {}
                    if(lf not in countAlign[j][k][le]):
                        countAlign[j][k][le][lf] = {}
                    countAlign[j][k][le][lf] = 0
                if(j not in totalAlign):
                    totalAlign[j] = {}
                if(le not in totalAlign):
                    totalAlign[j][le] = {}
                if(lf not in totalAlign):
                    totalAlign[j][le][lf] = 0


        totalF = {} #1 things, line 6 pseudo
        countEF = {} #2 things, line 5 pseudo
            #initializations
        for i in range(fileLength):
            for f in list(forSent[i].split()):
                totalF[f] = 0
                for e in list(engSent[i].split()):
                    if(e not in countEF):
                        countEF[e] = {}
                    countEF[e][f] = 0
    
        # line 9 of pseudocode
        tot = {}
        for i in range(fileLength):
            le = len(engSent[i].split())
            lf = len(forSent[i].split())


            # line 11 pseudocode -- compute normalization
            for j in range(1, len(engSent[i].split())):
                e = engSent[i].split()[j]
                #print(e)
                tot[e] = 0
                for f in list(forSent[i].split()):
                    #print(f)
                    tot[e] = pfe[e][f] * alignModel2[f][e][le][lf]


            # line 18 pseudocode -- collect counts 
            for e in range(1, len(engSent[i].split())):
                j = engSent[i].split()[e]
                for k in list(forSent[i].split()):
                    deltaThing = pfe[j][k] * alignModel2[k][j][le][lf] / tot[j]
                    countEF[j][k] += deltaThing
                    totalF[k] += deltaThing
                    countAlign[k][j][le][lf] += deltaThing
                    totalAlign[k][le][lf] += deltaThing


            #line 29 pseudo -- estimate probabilites


        # t(e|f) = 0 for all e, f -- line 30
        for e in engVocab: # collected "vocab" when reading in the sentences
            for f in forVocab:
                pfe[e][f]= 0
        
        # a(i| j, le, lf) = 0 for all i, j, le, lf -- line 31
        for i in alignModel2:
            for j in alignModel2[i]:
                for le in alignModel2[i][j]:
                    for lf in alignModel2[i][j][le]:
                        alignModel2[i][j][le][lf] = 0


        # for all e, f do t(e|f) = count(e|f)/total(f) -- line 32 & 33
        for e in engVocab:
            for f in forVocab:
                if e in pfe and f in pfe[e] and e in countEF and f in countEF[e] and f in totalF:
                    pfe[e][f] = countEF[e][f]/totalF[f]

        # line 35 
        for i in range(0, len(engSent)):
            lf = len(forSent[i].split())
            le = len(engSent[i].split())
            for j in range(0, lf):
                for k in range(1, le):
                    f = forSent[i].split()[j]
                    e = engSent[i].split()[k]
                    alignModel2[f][e][le][lf] = countAlign[f][e][le][lf] / totalAlign[f][le][lf]

        
     print(alignModel2)
     print("pfe")
     print(pfe)



     print("total Align")
   #   print(totalAlign)
     # for all i, j, le, lf do a(i| j, le, lf) = countA(i|j, le, lf)/totalA(i|j, le, lf)
     for i in alignModel2:
       for j in alignModel2[i]:
           for le in alignModel2[i][j]:
               for lf in alignModel2[i][j][le]:
                   county = countAlign[i][j][le][lf]
                   aligny = totalAlign[i][le][lf]
                   alignModel2[i][j][le][lf] = county/aligny


     # print p(f|e) table sorted by english and foreign words
     # print("\n"+e_sent+" "+f_sent+" p(f|e) Table")
  
     if fileLength > 10:
         rand10 = {}
         for _ in range(30):
             eng, forChoices = random.choice(list(pfe.items()))
             rand10[eng] = forChoices
         for eng, forWords in sorted(rand10.items(), key=lambda x: x[0]):
             for key, val in sorted(forWords.items(), key=lambda x: x[0]):
                 if val > threshhold:
                     print(eng+"\t"+key+"\t"+str(val))
     else:
         for eng, forWords in sorted(pfe.items(), key=lambda x: x[0]):
             for key, val in sorted(forWords.items(), key=lambda x: x[0]):
                 highestProb = 0
                 if val > threshhold:
                     print(eng+"\t"+key+"\t"+str(val))
      


"""
     # generate alignments if needed
     if(alignmentsflag == True):
         picks = []
         if(fileLength > 10):
             # pick 10 random sentences
             for _ in range(10):
                 picks.append(random.randint(0, fileLength-1))
         else:
             for i in range(fileLength):
                 picks.append(i)
      
         for idx in picks:
             print("\n"+"Foreign Sentence: "+forSent[idx]+"\n"+"English Sentence: "+engSent[idx])
             # for each foreign word
             for foreign in forSent[idx].split():
                 # go through all english words in sentence
                 ForToEngidx = -1
                 highestProb = 0
                 engWords = engSent[idx].split()
                 ForToEngWord = "UNK"
                 for j in range(len(engWords)):
                     # pick max of all the p(f|e) vals
                     if(pfe.get(engWords[j])[foreign] > highestProb):
                         highestProb = pfe.get(engWords[j])[foreign]
                         ForToEngidx = j
                         ForToEngWord = engWords[j]
                 print(foreign+" -> "+ForToEngWord)
"""


# requires 5th arg (boolean) stating if you want alignments for this corpus
#test0 = em('test.en 1', 'test.es 1', 20, 0.1, True)
#test1 = em('test.en', 'test.fr', 10, 0, True)
#fullRun = em('es-en.10k.en', 'es-en.10k.es', 10, 0.3, True)
test2 = ibm2("test.en", "test.fr", 4, 0.3, True)


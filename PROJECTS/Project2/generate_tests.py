#!/usr/bin/env python
# coding: utf-8

# In[1]:


#this script makes ten test cases for the citations and meta txt files inorder to quickly test parsing
#todo add argument parsing for meta and citations and error proofing
import os


# In[2]:



meta=open("arxiv-metadata.txtDGFP","r").readlines()
citations=open("arxiv-citations.txt","r").readlines()


# In[3]:


cite_dict={}
for i in range(len(citations)):
    citations[i]=citations[i].replace('\n','')
    if citations[i]=="+++++":
        k=i-1
        cites=[]
        while citations[k]!="-----":
            cites.append(citations[k])
            k-=1
        paper=citations[k-1]
        cite_dict[paper]=cites


# In[4]:


citdir='citation_tests/'
if not os.path.exists(citdir):
    os.mkdir(citdir)
for i,paper in enumerate(cite_dict.keys()):
    if i%10000==0:
        if i==100000:
            break
        testfile=open(citdir+'citation_test_{}.txt'.format(int(i/10000)),'a+')
    testfile.write(paper+'\n')
    testfile.write('-----\n')
    for link in cite_dict[paper]:
        testfile.write(link+'\n')
    testfile.write('+++++\n')
    i+=1


# In[5]:


metadir='meta_tests/'
if not os.path.exists(metadir):
    os.mkdir(metadir)
for i in range(0,len(meta),5):
    if i%10000==0:
        if i==100000:
            break
        testfile.close()
        testfile=open(metadir+'meta_test_{}.txt'.format(int(i/10000)),'a+')
    paper=meta[i:i+5]
    paper=str.join('',paper)  
    testfile.write(paper)
    


# In[ ]:





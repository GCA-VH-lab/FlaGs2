import argparse
import random
from random import randint
import colorsys
import os, sys, os.path, math
from tkinter import *
import csv
master = Tk()

usage= ''' Description:  Generates FlaGs tkinter based graphics; Requirement= Python3, BioPython; tkinter ; Optional Requirement= ETE3. '''
parser = argparse.ArgumentParser(description=usage)

parser.add_argument("-i", "--input", help=" TreeOrder_operon FlaGs output eg. GCF_accesion_output_TreeOrder_operon.tsv. ")
parser.add_argument("-ts", "--tshape", help=" Size of triangle shapes that represent flanking genes. Default = 12 ")
parser.add_argument("-tf", "--tfontsize", help=" Size of font inside triangles that represent flanking genes. Default = 4 ")
parser.add_argument("-o", "--out_prefix", required= True, help=" Any Keyword to define your output eg. MyQuery. Output file would be MyQuery_flankgenes.ps ")
args = parser.parse_args()
parser.parse_args()

from tkinter.font import Font #Font for postscript to-scale output
myFont12 = Font(family="Helvetica", size=12)
myFont7 = Font(family="Helvetica", size=7)


def postscriptSize(item):
    if int(item)<1000:
        return(0)
    else:
        return(int(item)/1000)
newQ=0

def operonFamily(item):
	if item==0:
		return ' '
	elif item==center:
		return ' '
	elif item==noProt:
		return ' '
	elif item==noProtP:
		return ' '
	elif item==noColor:
		return ' '
	else:
		return item

#Color generator
def random_color(h=None):
	"""Generates a random color in RGB format."""
	if not h:
		c = random.random()
	d = 0.5
	e = 0.5
	return _hls2hex(c, d, e)

def _hls2hex(c, d, e):
	return '#%02x%02x%02x' %tuple(map(lambda f: int(f*255),colorsys.hls_to_rgb(c, d, e)))

def outliner (item):
	if item =='#ffffff':
		return '#bebebe'
	elif item =='#f2f2f2':
		return '#008000'
	elif item =='#f2f2f3':
		return '#000080'
	else:
		return item

color={}
colorDict={}

if args.tshape:
	if int(args.tshape)>0:
		size=int(args.tshape)
	else:
		print("Please input the size of triangles, recommended 12. Not applicable for 0 and negative values")
else:
	size=12

if args.tfontsize:
	if int(args.tfontsize)>0:
		fsize=str(args.tfontsize)
	else:
		print("Please input the font Size required inside triangles, recommended 4. Not applicable for 0 and negative values")
else:
	fsize=str(4)

myFont12 = Font(family="Helvetica", size=size)
myFont7 = Font(family="Helvetica", size=fsize)

nPos=[]
pPos=[]


# tsvfile=open(args.input,'r', encoding='utf8').read()
with open(args.input, 'r', encoding='utf8') as tsvfile:
    tsvreader = csv.reader(tsvfile, delimiter="\t")
    qID=""
    for line in tsvreader:
        if len(line)>2:
            nP=line[5]
            pP=line[6]
            nPos.append(nP)
            pPos.append(pP)
            if qID != line[0]:
                newQ+=1
                qID=line[0]
                


eg1=open(args.input,'r').read()

windowMost=round(((int(max(pPos))+abs(int(min(nPos)))+1)*4)/100)
widthM=(windowMost*10)+750
heightM=int(newQ)*20

canvas = Canvas(master, width=widthM,height=heightM,background='white', scrollregion=(0,0,round(widthM*2.5),round(heightM*2.5)))
hbar=Scrollbar(master,orient=HORIZONTAL)
hbar.pack(side=BOTTOM,fill=X)
hbar.config(command=canvas.xview)
vbar=Scrollbar(master,orient=VERTICAL)
vbar.pack(side=RIGHT,fill=Y)
vbar.config(command=canvas.yview)
canvas.config(xscrollcommand=hbar.set, yscrollcommand=vbar.set)
canvas.pack(side=LEFT,expand=True,fill=BOTH)



egs=eg1.split("\n\n\n\n")
line_pos_y=0
for eg in egs:
    if eg!='':
        coln=0
        entries=eg.splitlines()
        ndoms=len(entries)
        ptnstats=entries[0].split("\t") #c2
        org=ptnstats[0][:ptnstats[0].index('|')]+ptnstats[0][ptnstats[0].index('|'):].replace('_',' ')
        textspace=widthM/2
        line_pos_y=line_pos_y + 16 - round(postscriptSize(newQ))
        half_dom_height=5-round(postscriptSize(newQ))
        text = canvas.create_text(textspace/2-textspace/8,line_pos_y, text=org, fill="#404040", font=myFont12)
        for entry in entries:
            items=entry.split("\t")
            aln_start=round(int(items[5])*4/100)
            aln_end=round(int(items[6])*4/100)
            strandType=items[3]
            dom1_name=int(items[4])
            dom1_len=(aln_end-aln_start)
            oL80=round(dom1_len*80/100)
            dom1_start=aln_start+textspace
            dom1_end=dom1_len+dom1_start
            id1=str(items[9])
            gene_start=int(items[5])


            # Colours (imported from FlaGs script)
            center=int(dom1_name)+1
            noProt=int(dom1_name)+2
            noProtP=int(dom1_name)+3
            noColor=int(dom1_name)+4

            color[noColor]='#ffffff'
            color[center]='#000000'
            color[noProt]='#f2f2f2'
            color[noProtP]='#f2f2f3'


            if dom1_name == 0:
                colorDict[dom1_name]=str('#ffffff')
            elif gene_start == 1:                    
                colorDict[dom1_name]=str('#000000')
            elif 'pseudogene_' in id1:
                colorDict[dom1_name]=str('#f2f2f2') 
            elif 'tRNA_' in id1:
                colorDict[dom1_name]=str('#f2f2f3')
            else:
                if dom1_name not in colorDict:
                    colorDict[dom1_name] = random_color()

            if strandType=='+':
                rect = canvas.create_polygon(dom1_start, line_pos_y+half_dom_height, dom1_start, line_pos_y-half_dom_height,dom1_start+oL80, line_pos_y-half_dom_height, dom1_end, line_pos_y, dom1_start+oL80, line_pos_y+half_dom_height, fill=colorDict[dom1_name], outline=outliner(colorDict[dom1_name]))
            else:
                rect = canvas.create_polygon(dom1_end-oL80, line_pos_y+half_dom_height, dom1_start, line_pos_y, dom1_end-oL80, line_pos_y-half_dom_height,dom1_end, line_pos_y-half_dom_height, dom1_end, line_pos_y+half_dom_height, fill=colorDict[dom1_name], outline=outliner(colorDict[dom1_name]))
            textd1 = canvas.create_text(dom1_start+(dom1_len/2),line_pos_y, text=operonFamily(dom1_name), font=myFont7, fill='black')
            coln=coln+1

retval = canvas.postscript(file=args.out_prefix+"_flankgenes.ps", height=heightM, width=widthM, colormode="color")



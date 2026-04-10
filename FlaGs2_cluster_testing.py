__author__		= "Fedor Taratorkin, Chayan Kumar Saha, Jose Nakamoto, Gemma C. Atkinson"
__copyright__	= "GNU General Public License v3.0"
__email__		= "taratorkinfed@gmail.com, jose.n.kuahara@gmail.com, chayan.sust7@gmail.com"

import re
import math
import argparse
import ftplib
import socket
import random
import time
import colorsys
import os, sys
import gzip
from collections import OrderedDict, Counter, defaultdict
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

import Bio
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import matplotlib as mpl

from collections import Counter  
import ete3

os.environ['QT_QPA_PLATFORM']='offscreen'

def parse_arguments():
    usage= ''' Description:  Identify flanking genes and cluster them based on similarity and visualize the structure; Requirement= Python3, BioPython. '''
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument("-a", "--assemblyList", help=" Protein Accession with assembly Identifier eg. GCF_000001765.3 in a text file separated by newline. ")
    parser.add_argument("-p", "--proteinList", help=" Protein Accession eg. XP_ or WP_047256880.1 in a text file separated by newline. ")
    parser.add_argument("-l", "--localGenomeList", help=" Genome File name and Protein Accession ")
    parser.add_argument("-ld", "--localGenomeDirectory", help=" Path for Local Files, Default directory is './' which is the same directory where the script is located or running from. ")
    parser.add_argument("-r", "--redundant", help=" To search all assembly type -r A or -r a but for selected number of assembly eg.,5 for each query use -r 5. ")
    parser.add_argument("-e", "--ethreshold", help=" E value threshold. Default = 1e-10 ")
    parser.add_argument("-n", "--number", help=" Number of Jackhmmer iterations. Default = 3")
    parser.add_argument("-g", "--gene", help=" Number of genes for looking up or downstream. Default = 4 ")
    parser.add_argument("-t", "--tree", action="store_true", help=" If you want to see flanking genes along with phylogenetic tree, requires ETE3 installation. By default it will not produce. ")
    parser.add_argument("-ts", "--tshape", help=" Size of triangle shapes that represent flanking genes, this option only works when -t is used. Default = 12 ")
    parser.add_argument("-tf", "--tfontsize", help=" Size of font inside triangles that represent flanking genes, this option only works when -t is used. Default = 4 ")
    parser.add_argument("-to", "--tree_order", action="store_true", help=" Generate Output with Tree, and then use the tree order to generate other view. ")
    parser.add_argument("-u", "--user_email", required=True, action="append", metavar="RECIPIENT",default=[], dest="recipients", help=" User Email Address (at least one required) ")
    parser.add_argument("-api", "--api_key", help="NCBI API Key, To get this key kindly check https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/ ")
    parser.add_argument("-o", "--out_prefix", required= True, help=" Any Keyword to define your output eg. MyQuery ")
    parser.add_argument("-c", "--cpu", help="Maximum number of parallel CPU workers to use for multithreads. ")
    parser.add_argument("-db", "--hmmdb", help=" Input hmm database to enable domain search in queries and flanking genes through hmmscan. Eg. Pfam-A.hmm, cd_all.hmm")
    #adding cutoff threshhold
    parser.add_argument("-k", "--keep", action="store_true", help=" If you want to keep the intermediate files eg. gff3 use [-k]. By default it will remove. ")
    parser.add_argument("-v", "--version", action="version", version='%(prog)s 1.1.5')
    parser.add_argument("-vb", "--verbose", action="store_true", help=" Use this option to see the work progress for each query as stdout. ")
    args = parser.parse_args()
    return args

def check_arguments(args):
	if args.proteinList and args.assemblyList:
		print ("Argument invalid: -a and -p can't be used together")
		sys.exit()
	if args.redundant:
		if not (int(args.redundant)>0 or args.redundant.lower()=='a'):
			print ("Argument invalid: -r value is not correct")
			sys.exit()
		if not args.proteinList:
			print ("Argument invalid: -r require -p")
			sys.exit()
	if args.cpu:
		if not int(args.cpu)>0:
			print ("Argument invalid: -c value is not correct")
			sys.exit()
	if not args.tree:
		if args.tshape or args.tfontsize or args.tree_order:
			print ("Argument invalid: -ts, -tf and -to require -t")
	else:
		if args.tshape:
			if not int(args.tshape)>0:
				print("Argument invalid: -ts value is not correct")
				sys.exit()
		if args.tfontsize:
			if not int(args.tfontsize)>0:
				print("Argument invalid: -tf value is not correct")
				sys.exit()
	if args.gene:
		if not int(args.gene)>0:
			print("Argument invalid: -g value is not correct")
			sys.exit()
	if args.localGenomeList:
		if args.localGenomeDirectory:
			if not os.path.isdir(args.localGenomeDirectory):
				print('No directory Found as : '+ args.localGenomeDirectory)
				sys.exit()
	else:
		if args.localGenomeDirectory:
			print("Arguments invalid: Please use -l flag to make -ld flag working")
			sys.exit()

#small utility functions

def checkBioPython(): #Checking Biopython Version
	return Bio.__version__
def random_color(h=None):
	if not h:
		c = int((random.randrange(0,100,5))*3.6)/100
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

def checkChar(item): #removing characters
	items=item.replace('\t','').replace(' ','')
	return re.sub("[a-zA-Z0-9_.]","",items)

def remBadChar(item): #removing characters from species name
	return re.sub("[^a-zA-Z0-9]"," ",item).replace(" ","_")

def des_check(item):
	if item:
		return item
	else:
		return 'notFound'

def normalize_strand(item1, item2):  #Strand direction change
	if item1=='+':
		return item2
	else:
		if item2=='+':
			return '-'
		else:
			return '+'

def up(item):
	if item=='+':
		return 'Upstream '
	else:
		return 'Downstream '

def down(item):
	if item=='+':
		return 'Downstream '
	else:
		return 'Upstream '

def ups(item):
	if item=='+':
		return '-'
	else:
		return '+'

def downs(item):
	if item=='+':
		return '+'
	else:
		return '-'

def lcheck(item):
	if 1 in item:
		return 1
	else:
		return 0

def postscriptSize(item):
	if int(item)<1000:
		return(0)
	else:
		return(int(item)/1000)


def getSpeciesFromGCF(faa,item):
	if item!='':
		if args.redundant:
			return remBadChar(item)+'_'+remBadChar(faa)
		else:
			return remBadChar(item)
	else:
		return 'Nothing'
def getGeneId(item):
	matchObject = re.search('(GeneID:.*?,)', item)
	if matchObject:
		return matchObject.group(1)[:-1]

def getGeneId_gene(item):
	matchObject2 = re.search('(GeneID:.*;)', item)
	if matchObject2:
		return matchObject2.group(1).split(';')[0]

def lenChecker(List):
	if List:
		found=0
		for item in List:
			if item==True and item>5000:
				found+=1
		if found==0:
			return 'keep'
		else:
			return 'delete'
	else:
		return 'keep'
def similarityID (item1, item2):
	if item1==item2:
		return 'Same'
	else:
		return 'Changed'

#spLocalrelated

def spLocal(faa,acc): #getting species name using assembly number or accession
	if faa in speciesNameFromOnlineDict:
		return speciesNameFromOnlineDict[faa]
	else:
		faaFile=faa+'.faa.gz'
		fastaSeq = gzip.open(localDir+faaFile, "rt")
		for record in SeqIO.parse(fastaSeq, "fasta"):
			if record.id==acc:
				if args.redundant:
					return remBadChar(record.description.split('[')[-1][:-1])+'_'+remBadChar(faa)
				else:
					return remBadChar(record.description.split('[')[-1][:-1])

def desLocal(faa,acc):
	faaFile=faa+'.faa.gz'
	fastaSeq = gzip.open(localDir+faaFile, "rt")
	for record in SeqIO.parse(fastaSeq, "fasta"):
		if record.id==acc:
			return record.description.split('[')[0]

def seqLocal(faa,acc):
	faaFile=faa+'.faa.gz'
	fastaSeq = gzip.open(localDir+faaFile, "rt")
	for record in SeqIO.parse(fastaSeq, "fasta"):
		if record.id==acc:
			return str(record.seq)

def localNone(item):
	if item==None:
		return '--'
	else:
		return item

def seqFasLocal(faa,acc): #making fasta file from accession
	if spLocal(faa,acc):
		faaFile=faa+'.faa.gz'
		fastaSeq = gzip.open(localDir+faaFile, "rt")
		for record in SeqIO.parse(fastaSeq, "fasta"):
			if record.id==acc.split('#')[0]:
				record.id=acc+'|'+spLocal(faa,acc)#.replace(':','_').replace('[','_').replace(']','_')
				record.description=''
				return record.format("fasta")
	else:
		faaFile=faa+'.faa.gz'
		fastaSeq = gzip.open(localDir+faaFile, "rt")
		for record in SeqIO.parse(fastaSeq, "fasta"):
			if record.id==acc.split('#')[0]:
				if args.redundant:
					record.id=acc+'|'+remBadChar(record.description.split('[')[-1][:-1])+'_'+remBadChar(faa)#.replace(':','_').replace('[','_').replace(']','_')
				else:
					record.id=acc+'|'+remBadChar(record.description.split('[')[-1][:-1])
				record.description=''
				return record.format("fasta")

def seqFasLenLocal(faa,acc): #making fasta file from accession
	if spLocal(faa,acc):
		faaFile=faa+'.faa.gz'
		fastaSeq = gzip.open(localDir+faaFile, "rt")
		for record in SeqIO.parse(fastaSeq, "fasta"):
			if record.id==acc.split('#')[0]:
				record.id=acc+'|'+spLocal(faa,acc)#.replace(':','_').replace('[','_').replace(']','_')
				record.description=''
				if record.seq:
					return len(record.seq)

	else:
		faaFile=faa+'.faa.gz'
		fastaSeq = gzip.open(localDir+faaFile, "rt")
		for record in SeqIO.parse(fastaSeq, "fasta"):
			if record.id==acc.split('#')[0]:
				if args.redundant:
					record.id=acc+'|'+remBadChar(record.description.split('[')[-1][:-1])+'_'+remBadChar(faa)#.replace(':','_').replace('[','_').replace(']','_')
				else:
					record.id=acc+'|'+remBadChar(record.description.split('[')[-1][:-1])
				record.description=''
				if record.seq:
					return len(record.seq)

def redundantCreate(setDict,nums):
	if nums=='A' or nums=='a':
		newList=random.sample(setDict,len(setDict))
	else:
		if len(setDict)>int(nums):
			newList=random.sample(setDict,int(nums))
		else:
			newList=random.sample(setDict,len(setDict))
	return newList

#Download assembly summary report from NCBI Refseq and genBank
def database_loader():
	refDb='./refSeq.db'
	genDb='./genBank.db'
	refDbSize, genDbSize = '0','0'
	if os.path.isfile(refDb):
		refDbSize=os.path.getsize(refDb)
	if os.path.isfile(genDb):
		genDbSize=os.path.getsize(genDb)

	ftp = ftplib.FTP('ftp.ncbi.nlm.nih.gov', 'anonymous', 'anonymous@ftp.ncbi.nih.gov')
	ftp.cwd("/genomes/refseq") # move to refseq directory

	filenames = ftp.nlst() # get file/directory names within the directory
	if 'assembly_summary_refseq.txt' in filenames:
		ftp.sendcmd("TYPE i")
		if int(ftp.size('assembly_summary_refseq.txt'))!=int(refDbSize):#check if the previously downloaded db exists and if that's updated to recent one
			ftp.retrbinary('RETR ' + 'assembly_summary_refseq.txt', open('refSeq.db', 'wb').write) # get the assembly summary from refseq

	ftp_gen = ftplib.FTP('ftp.ncbi.nlm.nih.gov', 'anonymous', 'anonymous@ftp.ncbi.nih.gov')
	ftp_gen.cwd("/genomes/genbank") # move to genbank directory

	filenames = ftp_gen.nlst() # get file/directory names within the directory
	if 'assembly_summary_genbank.txt' in filenames:
		ftp_gen.sendcmd("TYPE i")
		if int(ftp_gen.size('assembly_summary_genbank.txt'))!=int(genDbSize):#check if the previously downloaded db exists and if that's updated to recent one
			ftp_gen.retrbinary('RETR ' + 'assembly_summary_genbank.txt', open('genBank.db', 'wb').write) # get the assembly summary from refseq

def query_list_builder(queryFile, expect_single_column):
    queryList=[]
    with open (queryFile, 'r') as qList:
        for line in qList:
            if checkChar(line.rstrip().replace(' ',''))=='':
                Line=line.rstrip().replace(' ','').split('\t')
                if expect_single_column and len(Line)==1:
                    queryList.append(Line)
                elif len(Line)>1 and not expect_single_column:
                    newFormat=Line[1]+'\t'+Line[0]
                    queryList.append(newFormat.split('\t'))
                else:
                    print('Check Input file, Incorrect Format.')
                    sys.exit()
            else:
                print('The submitted query might include characters not found in NCBI protein accessions eg. > , # , ! etc. Please provide correct format, Thanks!')
                sys.exit()
    return queryList

#entrez communication 

def batch_accession_from_wp(accession_list):
	accession_list = list(accession_list)
	if not accession_list:
		return {}
	results = {acc: False for acc in accession_list}
	i = 1
	retry = True
	while retry and i < 6:
		try:
			variableIsquare = i**2
			changedtimeout = 10*variableIsquare
			socket.setdefaulttimeout(changedtimeout)
			time.sleep(ncbi_time)
			handle = Entrez.efetch(db="ipg", id=','.join(accession_list),
								   rettype="ipg", retmode="text")
			if handle:
				retry = False
				found = {acc: set() for acc in accession_list}
				for item in handle:
					if item[0:2] != 'Id':
						if re.search("GC._\d*\.\d", item):
							itemLine = item.rstrip().split('\t')
							iAccession = itemLine[6]
							iAssembly  = itemLine[-1]
							if iAccession in found:
								found[iAccession].add(iAssembly)
				handle.close()
				for acc in accession_list:
					results[acc] = found[acc] if found[acc] else {'NAI'}
			else:
				i += 1
		except Exception as e:
			retry = True
			i += 1
			if not i < 6:
				print("\t\tQuery {}, not found in database. \n"
					  "\t\tContinuing with the next protein in the list ... \n".format(
					  ', '.join(accession_list)))
	return results

def sortGCFvsGCA(gcagcfSet):
	if gcagcfSet!='NAI':
		Aset=set()
		Fset=set()
		for items in gcagcfSet:
			if items[2]=='A':
				Aset.add(items)
			if items[2]=='F':
				Fset.add(items)
		if len(Fset)>0:
			return Fset
		elif len(Fset)==0 and len(Aset)>0:
			return Aset
		else:
			return gcagcfSet
	else:
		return gcagcfSet

def batch_accession_from_xp(accession_list):
	accession_list = list(accession_list)
	if not accession_list:
		return {}
	results = {acc: False for acc in accession_list}
	i = 1
	retry = True
	while retry and i < 6:
		try:
			variableIsquare = i**2
			changedtimeout = 10*variableIsquare
			socket.setdefaulttimeout(changedtimeout)
			time.sleep(ncbi_time)
			handle = Entrez.efetch(db="protein", id=','.join(accession_list),
								   rettype="gbwithparts", retmode="text")
			if handle:
				retry = False
				records = list(SeqIO.parse(handle, "genbank"))
				handle.close()
				for record in records:
					rec_id   = record.id
					rec_base = rec_id.split('.')[0]
					matched  = None
					for acc in accession_list:
						if acc == rec_id or acc == rec_base:
							matched = acc
							break
					if matched is None:
						continue
					bio = {ref.split(':')[1]
						   for ref in record.dbxrefs
						   if ref.split(':')[0] == 'BioProject'}
					results[matched] = bio if bio else {'NAI'}
			else:
				i += 1
		except Exception as e:
			retry = True
			i += 1
			if not i < 6:
				print("\t\tQuery {}, not found in database. \n"
					  "\t\tContinuing with the next protein in the list ... \n".format(
					  ', '.join(accession_list)))
	return results

def accession_from_xp(accession_nr):
	return batch_accession_from_xp([accession_nr]).get(accession_nr, False)

def accession_from_wp(accession_nr):
	return batch_accession_from_wp([accession_nr]).get(accession_nr, False)

def seq_from_wp(accession_nr):
	if accession_nr[-1]!='*':
		#Entrez.email = "_@gmail.com"  # If you do >3 entrez searches on NCBI per second, your ip will be
		# blocked, warning is sent to this email.
		try:
			time.sleep(ncbi_time)
			handle = Entrez.efetch(db="protein", id=accession_nr, rettype="gbwithparts", retmode="text")
			
		except Exception as e:
			print(str(e), ", error in entrez-fetch protein accession, {}, not found in database. \n" "Continuing with the next protein in the list. \nError in function: {}".format(accession_nr, seq_from_wp.__name__))
			return False

		record = SeqIO.read(handle, "genbank")
		handle.close()
		return record.description.split('[')[0]+'\t'+record.seq
	else:
		return accession_nr[:-1]+'\t'+'--'

def _fetch_ipg_text(accnr):
	i = 1
	while i < 6:
		try:
			socket.setdefaulttimeout(10 * i**2)
			time.sleep(ncbi_time)
			epost_1 = Entrez.read(Entrez.epost(db="protein", id=accnr))
			return Entrez.efetch(db="protein", rettype='ipg', retmode='text',
								 webenv=epost_1["WebEnv"],
								 query_key=epost_1["QueryKey"])
		except Exception:
			i += 1
	raise RuntimeError("Failed to fetch IPG for {}".format(accnr))

def identicalProtID(accnr): #searching for identical proteins
	try:
		iden_prots = _fetch_ipg_text(accnr)
	except RuntimeError:
		return accnr
	iAccSet=set()
	sAccSet=set()
	inrAccSet=set()
	snrAccSet=set()
	for item in iden_prots:
		if item[0:2]!='Id':
			if re.search("GC._\d*\.\d", item):
				itemLine=item.rstrip().split('\t')
				iAccession=itemLine[6]
				if iAccession!=accnr and iAccession[2]=='_':
					iAccSet.add(iAccession)
				if iAccession!=accnr and iAccession[2]!='_':
					inrAccSet.add(iAccession)
				if iAccession==accnr and iAccession[2]=='_':
					sAccSet.add(iAccession)
				if iAccession==accnr and iAccession[2]!='_':
					snrAccSet.add(iAccession)
	if len(iAccSet)>0 and len(sAccSet)==0:
		return random.sample(list(iAccSet),1)[0]
	elif len(sAccSet)!=0:
		return random.sample(list(sAccSet),1)[0]
	elif len(inrAccSet)>0 and len(snrAccSet)==0 and len(iAccSet)==0:
		return random.sample(list(inrAccSet),1)[0]
	else:
		return accnr

def identicalProtID_WP(accnr): #searching for identical proteins
	try:
		iden_prots = _fetch_ipg_text(accnr)
	except RuntimeError:
		return accnr
	iAccSet=set()
	for item in iden_prots:
		if item[0:2]!='Id':
			if re.search("GC._\d*\.\d", item):
				itemLine=item.rstrip().split('\t')
				iAccession=itemLine[6]
				if iAccession!=accnr and iAccession[:3]=='WP_':
					iAccSet.add(iAccession)
	if len(iAccSet)>0:
		return random.sample(list(iAccSet),1)[0]
	else:
		return accnr


def identicalProtID_WP_Sp(accnr): #searching for identical proteins with same assembly  for 'NP_417570.1' > WP_000785722.1|GCF_000005845.2
	try:
		iden_prots = _fetch_ipg_text(accnr)
	except RuntimeError:
		return '#'
	iAccSetSpecial=set()
	iAccNRSetSpecial=set()
	iAssemblyList=[]
	iAssemblyListNR=[]
	for item in iden_prots:
		if item[0:2]!='Id':
			if re.search("GC._\d*\.\d", item):
				itemLine=item.rstrip().split('\t')
				iAccession=itemLine[6]
				iAssembly=itemLine[-1]
				if iAccession==accnr and iAccession[2]=='_':
					iAssemblyList.append(iAssembly)
				if iAccession==accnr and iAccession[-2]=='.' and iAssembly[0]=='G':
					iAssemblyListNR.append(iAssembly)
				if iAccession!=accnr and iAccession[-2]=='.' and iAssembly[0]=='G':
					iAccNRSetSpecial.add(iAccession+'|'+iAssembly)
				if iAccession!=accnr and iAccession[:3]=='WP_':
					if iAssemblyList:
						if iAssembly==iAssemblyList[0]:
							iAccSetSpecial.add(iAccession+'|'+iAssembly)
	if iAccSetSpecial:
		return random.sample(list(iAccSetSpecial),1)[0]
	elif iAssemblyList:
		return accnr+'|'+random.sample(iAssemblyList,1)[0]
	elif iAssemblyListNR:
		return accnr+'|'+random.sample(iAssemblyListNR,1)[0]
	elif iAccNRSetSpecial:
		return random.sample(list(iAccNRSetSpecial),1)[0]
	else:
		return '#'

def identicalProtID_redundant(accnr): #searching for identical proteins with same assembly  for 'NP_417570.1' > WP_000785722.1|GCF_000005845.2
	try:
		iden_prots = _fetch_ipg_text(accnr)
	except RuntimeError:
		return '#'
	iAssemblyset=set()
	for item in iden_prots:
		if item[0:2]!='Id':
			if re.search("GC._\d*\.\d", item):
				itemLine=item.rstrip().split('\t')
				iAccession=itemLine[6]
				iAssembly=itemLine[-1]
				if iAccession==accnr:
					iAssemblyset.add(iAssembly)
	return iAssemblyset if iAssemblyset else '#'

def reporter(i1,i2,i3,i4,i5):
	if i2=='No' and i3=='No' and i4=='No' and i5=='No':
		return i1 + ' Failed :  No record for accession ' + i1
	if i2!='No' and i3=='Same' and i4!='No' and i5=='No':
		return i1 + ' is a valid NCBI protein accession but Discarded :  No Flanking Gene was found for ' + i1 + ' in Assembly ID '+ i4
	if i2!='No' and i3=='Same' and i4!='No' and i5!='No':
		return i1 + ' is valid as a NCBI protein accession and reported in Assembly ID '+ i4
	if i2!='No' and i3=='Changed' and i4!='No' and i5=='No':
		return i1 + ' is invalid NCBI protein accession therefore converted to identical RefSeq sequence with accession '+ i2 + ' but Discarded :  No Flanking Gene was found for ' + i2 + 'in Assembly ID '+ i4
	if i2!='No' and i3=='Changed' and i4!='No' and i5!='No':
		return i1 + ' is invalid NCBI protein accession therefore converted to identical RefSeq sequence with accession '+ i2 + ' which is reported in Assembly ID '+ i4

def _download_assembly(query, item, query_idx, total_queries):

	if args.localGenomeList or item not in accnr_list_dict:
		return item, None, False

	species_label = getSpeciesFromGCF(item, accnr_list_dict[item].split('\t')[0])

	for attempt in range(1, 6):
		AssemDown  = 0
		AssemFailed = 0
		try:
			ftpLine      = accnr_list_dict[item].split('\t')[1]
			ftp_path     = '/'.join(ftpLine.split('/')[3:])
			ftp = ftplib.FTP('ftp.ncbi.nlm.nih.gov', 'anonymous', 'anonymous@ftp.ncbi.nih.gov')
			ftp.set_pasv(True)
			ftp.cwd('/' + ftp_path)
			files = ftp.nlst()
			FileToDownload = [f for f in files
							  if '_genomic.gff.gz' in f or '_protein.faa.gz' in f]

			if len(FileToDownload) == 2:
				for elements in FileToDownload:
					if args.verbose:
						ftp.set_debuglevel(1)
					ftp.voidcmd('TYPE I')
					ftp.sendcmd("TYPE i")
					is_gff = '_genomic.gff.gz' in elements
					out_name = item + ('.gff.gz' if is_gff else '.faa.gz')
					try:
						remote_size = ftp.size(elements)
						with open(out_name, 'wb') as fout:
							ftp.retrbinary('RETR ' + elements, fout.write)
						ftp.sock.setsockopt(socket.SOL_SOCKET,   socket.SO_KEEPALIVE, 1)
						ftp.sock.setsockopt(socket.IPPROTO_TCP,  socket.TCP_KEEPINTVL, 15)
						ftp.sock.setsockopt(socket.IPPROTO_TCP,  socket.TCP_KEEPCNT,   8)
						ftp.voidcmd("NOOP")
						if os.path.isfile(localDir + out_name) and os.path.getsize(out_name) == remote_size:
							AssemDown += 1
					except Exception:
						AssemFailed -= 1
				ftp.close()
			else:
				AssemFailed -= 2
				ftp.close()
				break  # No point retrying if files aren't there

			if AssemDown + AssemFailed == 2:
				if args.verbose:
					with _print_lock:
						print('\n\t> NCBI Genome Assembly has been downloaded for '
							  + query.split('#')[0]
							  + ' (' + str(query_idx) + '/' + str(total_queries) + ')\n')
				return item, species_label, True

		except Exception:
			if attempt == 5:
				break
			continue  # retry

	return item, species_label, False

#output functions 

def write_operon_tsv(out_filename, query_order):
	"""Write operon TSV and return (nPos, pPos) coordinate lists."""
	nPos=[]
	pPos=[]
	with open(out_filename, 'w') as opOut:
		for queries in query_order:
			for items in sorted(accFlankDict[queries]):
				ids=accFlankDict[queries][items][:-1]
				lengths=LengthDict[ids]
				species=queries+'|'+remBadChar(speciesDict[queries])
				qStrand=queryStrand[queries]
				nStrand=accFlankDict[queries][items][-1]
				family=familyDict[ids.split('#')[0]]
				info=acc_CGF_Dict[queries] if queries in acc_CGF_Dict else 'not_found\tnot_found\tnot_found'
				if qStrand=='+':
					startPos=int(positionDict[accFlankDict[queries][0][:-1]].split('\t')[0])-1
					start=int(positionDict[ids].split('\t')[0])
					end=int(positionDict[ids].split('\t')[1])
					print(species, lengths, qStrand, nStrand, family, start-startPos, end-startPos, start, end, ids, info, sep='\t', file=opOut)
					nPos.append(start-startPos)
					pPos.append(end-startPos)
				else:
					startPos=int(positionDict[accFlankDict[queries][0][:-1]].split('\t')[1])+1
					start=int(positionDict[ids].split('\t')[1])
					end=int(positionDict[ids].split('\t')[0])
					print(species, lengths, qStrand, nStrand, family, startPos-start, startPos-end, end, start, ids, info, sep='\t', file=opOut)
					nPos.append(startPos-start)
					pPos.append(startPos-end)
			print('\n\n', file=opOut)
	return nPos, pPos

def draw_operon_pdf(tsv_filename, pdf_filename):
	"""Read an operon TSV and render it as a PDF."""
	mpl.rcParams['pdf.fonttype'] = 42
	mpl.rcParams['ps.fonttype'] = 42
	arrowList = []
	accession_List = []
	gene_start_list=[]
	gene_end_list=[]
	fig, ax = plt.subplots(1, 2, sharey='row', gridspec_kw={'width_ratios': [5, 15]})
	main_file = open(tsv_filename,'r').read()
	eg1 = main_file.split("\n\n\n\n")
	y_level_m = 0
	for m in eg1:
		if m != '':
			entries1 = m.splitlines()
			y_level_m = y_level_m-2.75
			for entry in entries1:
				items1 = entry.split("\t")
				x_gene_start = int(items1[5])
				x_gene_end = int(items1[6])
				dx_gene_length = int(items1[1])
				gene_direction = items1[3]
				dom1_name = int(items1[4])
				id1 = str(items1[9])
				accession = str(items1[0])
				accession_List.append(accession)
				gene_start_list.append(x_gene_start)
				gene_end_list.append(x_gene_end)
				# 3a. Keep arrow shape for very short genes
				if dx_gene_length < 100:
					dx_gene_length = int(items1[1])*2.5
				# 3b. Point arrow head in correct direction
				if gene_direction == '-':
					x_gene_start = int(items1[6])
					x_gene_end = int(items1[5])
					dx_gene_length = dx_gene_length*(-1)
				else:
					x_gene_start = int(items1[5])
					x_gene_end = int(items1[6])
				# 3c. Draw gene as arrow
				head_len = abs(dx_gene_length)/2 if abs(dx_gene_length) < 200 else 150
				arrowList.append(ax[1].arrow(x=x_gene_start, y=y_level_m*0.65, dx=dx_gene_length, dy=0,
					width=1.1, head_width=1.1, length_includes_head=True, head_length=head_len,
					facecolor=colorDict[dom1_name], edgecolor=outliner(colorDict[dom1_name]), alpha=1))
				# 5. Family number label inside arrow
				text_x = x_gene_start + (dx_gene_length/2)
				if dom1_name != 0 and 'pseudogene_' not in id1 and 'RNA_' not in id1 and 'other' not in id1 and x_gene_start != 1:
					ax[1].text(text_x, y_level_m*0.65-0.35, s=dom1_name, horizontalalignment='center', color='k', font={'family':'sans-serif','size':8})
				# 6. Organism label on left panel
				ptnstats = entries1[0].split("\t")
				org = ptnstats[0][:ptnstats[0].index('|')]+ptnstats[0][ptnstats[0].index('|'):].replace('_',' ')
				ax[0].text(0.5, y_level_m*0.65-0.35, org, horizontalalignment='center', color='#000000', font={'family':'sans-serif','size':10})
				ax[0].set_axis_off()
	y_new = (len(Counter(accession_List).keys()))*0.25
	x_new = (max(gene_end_list) - min(gene_start_list))/3000 + 10
	fig.set_figheight(y_new)
	fig.set_figwidth(x_new)
	plt.axis('off')
	plt.savefig(pdf_filename, format='pdf', bbox_inches='tight')
	plt.close()


starttime = time.perf_counter()
args=parse_arguments()
check_arguments(args)

#setting parameters
localDir='./'
if args.ethreshold:
	evthresh=args.ethreshold
else:
	evthresh="1e-10"
if args.number:
	iters=args.number
else:
	iters="3"
if args.tree:
	if args.tshape:
		size=args.tshape
	else:
		size=12
	if args.tfontsize:
		fsize=args.tfontsize
	else:
		fsize="4"
if args.gene:
	if int(args.gene)>0:
		s= str(int(args.gene)+1)
else:
	s=5
if not args.localGenomeList:
	if args.api_key:
		Entrez.api_key = args.api_key #Valid API-key allows 10 queries per seconds, which makes the tool run faster
else:
	if args.api_key:
		print('Since FlaGs2 will use Local Data api_key is not necessary, Thanks!')
		sys.exit()
core = int(args.cpu) if args.cpu else 1
print("\nStarting FlaGs2 version 1.1.5 \nPlease only run one instance of FlaGs2 at a time to avoid making more queries than NCBI’s limit.")
print('For more information, please check https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/ \n')
print('Checking for RefSeq and Genbank summary files and downloading if needed ... \n')

Entrez.tool = 'FlaGs2'
ncbi_time= 0.4
timeout = 10
socket.setdefaulttimeout(timeout)
Entrez.email = args.recipients[0] #User email
Entrez.max_tries = 5
Entrez.sleep_between_tries = 60

queryList=[]

if args.localGenomeList: queryList=query_list_builder(args.localGenomeList, False)
else:
	if args.proteinList: queryList=query_list_builder(args.proteinList, True)
	else:queryList=query_list_builder(args.assemblyList, False)
	database_loader()
	assemblyName={}
	bioDict={} #bioproject as keys and assemble number (eg.GCF_000001765.1) as value
	accnr_list_dict={} #create a dictionary accessionNumber is a key and Organism name and ftp Gff3 download Link as value
	with open('refSeq.db', 'r') as fileIn:
		for line in fileIn:
			if line[0]!='#':
				Line=line.rstrip().split('\t')
				accnr_list_dict[Line[0]]= Line[7]+'\t'+Line[19]
				bioDict[Line[1]]=Line[0]
				assemblyName[Line[0]]=Line[0]
	ftp_gen = ftplib.FTP('ftp.ncbi.nlm.nih.gov', 'anonymous', 'anonymous@ftp.ncbi.nih.gov')
	ftp_gen.cwd("/genomes/genbank") # move to refseq directory
	assemblyName_GCA={}
	bioDict_gen={}
	accnr_list_dict_gen={} #create a dictionary accessionNumber is a key and Organism name and ftp Gff3 download Link as value
	with open('genBank.db', 'r') as fileIn:
		for line in fileIn:
			if line[0]!='#':
				Line=line.rstrip().split('\t')
				if len(Line)>19:
					if Line[18]=='identical':
						if Line[17] in accnr_list_dict:
							bioDict_gen[Line[1]]=Line[0]
							accnr_list_dict_gen[Line[0]]= accnr_list_dict[Line[17]]
							assemblyName_GCA[Line[0]]=Line[17]
						else:
							accnr_list_dict_gen[Line[0]]=Line[7]+'\t'+Line[19]
					else:
						accnr_list_dict_gen[Line[0]]=Line[7]+'\t'+Line[19]

	bioDict.update(bioDict_gen)
	accnr_list_dict.update(accnr_list_dict_gen)
	assemblyName.update(assemblyName_GCA)
	ftp_gen.close()

	print ('\n'+ '>> Database Downloaded. Cross-checking of the accession list in progress ...'+ '\n')
q=0
ne=0
queryDict={} #protein Id as query and a set of assembly number as value [either All or Species of interest]
#{'WP_019504790.1#1': {'GCF_000332195.1'}, 'WP_028108719.1#2': {'GCF_000422645.1'}, 'WP_087820443.1#3': {'GCF_900185565.1'}}
#{'WP_019504790.1#1': 'GCF_000332195.1', 'WP_028108719.1#2': 'GCF_000422645.1', 'WP_087820443.1#3': 'GCF_900185565.1'} local
if not args.localGenomeList:
	with open (args.out_prefix+'_NameError.txt', 'w') as fbad:

		# Classify queries by type for batching
		single_indices       = [i for i, qr in enumerate(queryList) if len(qr) < 2]
		paired_indices       = [i for i, qr in enumerate(queryList) if len(qr) >= 2]

		wp_single_indices    = [i for i in single_indices
								if queryList[i][0][:2]=='WP' and queryList[i][0][-2]=='.']
		xp_single_indices    = [i for i in single_indices
								if queryList[i][0][:2]=='XP' and queryList[i][0][-2]=='.']
		other_single_indices = [i for i in single_indices
								if not (queryList[i][0][:2]=='WP' and queryList[i][0][-2]=='.')
								and not (queryList[i][0][:2]=='XP' and queryList[i][0][-2]=='.')]

		xp_paired_indices    = [i for i in paired_indices
								if queryList[i][0][:3]=='XP_' and queryList[i][0][-2]=='.']
		wp_paired_indices    = [i for i in paired_indices
								if queryList[i][0][:3]!='XP_' and queryList[i][0][-2]=='.']

		# Batch Entrez fetches 

		wp_single_results = (batch_accession_from_wp([queryList[i][0] for i in wp_single_indices])
							 if wp_single_indices else {})

		xp_single_results = (batch_accession_from_xp([queryList[i][0] for i in xp_single_indices])
							 if xp_single_indices else {})
		other_resolved = {}
		for i in other_single_indices:
			other_resolved[i] = identicalProtID(queryList[i][0])

		other_needs_wp = [i for i in other_single_indices
						  if other_resolved[i] != queryList[i][0]
						  and other_resolved[i][:-3] != 'XP_']
		other_needs_xp = [i for i in other_single_indices
						  if other_resolved[i] != queryList[i][0]
						  and other_resolved[i][:-3] == 'XP_']
		other_special  = [i for i in other_single_indices
						  if other_resolved[i] == queryList[i][0]]

		other_wp_results = {}
		if other_needs_wp:
			other_wp_results = batch_accession_from_wp(
				list({other_resolved[i] for i in other_needs_wp}))

		other_xp_results = {}
		if other_needs_xp:
			other_xp_results = batch_accession_from_xp(
				list({other_resolved[i] for i in other_needs_xp}))

		special_exceptional_wp = {}
		special_special_out    = {}
		for i in other_special:
			rid = other_resolved[i]
			special_exceptional_wp[i] = identicalProtID_WP(rid)
			special_special_out[i]    = identicalProtID_WP_Sp(rid)

		exceptional_wp_ids = list({v for v in special_exceptional_wp.values() if v[:3]=='WP_'})
		exceptional_wp_results = (batch_accession_from_wp(exceptional_wp_ids)
								  if exceptional_wp_ids else {})

		xp_paired_results = (batch_accession_from_xp([queryList[i][0] for i in xp_paired_indices])
							 if xp_paired_indices else {})
		wp_paired_results = (batch_accession_from_wp([queryList[i][0] for i in wp_paired_indices])
							 if wp_paired_indices else {})

		for idx, query in enumerate(queryList):
			q+=1
			accession_from_wp_out=''
			accession_from_xp_out=''
			identicalProtID_Out=''
			accession_from_wp_ID_out=''
			accession_from_xp_ID_out=''
			accession_from_wp_IDSame_out=''
			exceptionalWP_out=''
			accession_from_wp_exceptional=''
			special_out = ''
			assembly_from_identical=''
			if args.verbose:
				print('\t Checking Query '+ query[0] +' ....'+ '('+str(q)+'/'+str(len(queryList))+')')
			if len(query)<2:
				if query[0][:2]=='WP' and query[0][-2]=='.': #WP Accession full WP_000785722.1
					accession_from_wp_out=wp_single_results.get(query[0], False)
					if accession_from_wp_out:
						queryDict[query[0]+'#'+str(q)]=sortGCFvsGCA(accession_from_wp_out)
					else:
						ne+=1
						print(query[0], file= fbad)
				elif query[0][:2]=='XP' and query[0][-2]=='.': #XP Accession full XP_003256407.1
					accession_from_xp_out=xp_single_results.get(query[0], False)
					if accession_from_xp_out:
						assemList=[]
						for bioprojs in accession_from_xp_out:
							if bioprojs in bioDict:
								assemList.append(bioDict[bioprojs])
						if assemList:
							queryDict[query[0]+'#'+str(q)]=sortGCFvsGCA(set(assemList))
					else:
						ne+=1
						print(query[0], file= fbad)
				else: #other accessions can be XP_003256407 , WP_000785722 too
					identicalProtID_Out=other_resolved[idx]
					if identicalProtID_Out!=query[0]: #can be anything XP_ or WP_ or YP_ or NP_
						if identicalProtID_Out[:-3]!='XP_': #Not XPs
							accession_from_wp_ID_out=other_wp_results.get(identicalProtID_Out, False)
							if accession_from_wp_ID_out:
								asset=set()
								for elements in accession_from_wp_ID_out:
									asset.add(elements)
								if len(asset)>0:
									queryDict[identicalProtID_Out+'#'+str(q)+'.'+query[0]]=sortGCFvsGCA(asset)
							else:
								ne+=1
								print(query[0], file= fbad)
						if identicalProtID_Out[:-3]=='XP_': #if XPs
							accession_from_xp_ID_out=other_xp_results.get(identicalProtID_Out, False)
							if accession_from_xp_ID_out:
								assemList=[]
								for bioprojs in accession_from_xp_ID_out:
									if bioprojs in bioDict:
										assemList.append(bioDict[bioprojs])
								if assemList:
									queryDict[identicalProtID_Out+'#'+str(q)+'.'+query[0]]=sortGCFvsGCA(set(assemList))
							else:
								ne+=1
								print(query[0], file= fbad)
					elif identicalProtID_Out==query[0]: #excluding XP Wp pre  #YP NP
						exceptionalWP_out = special_exceptional_wp[idx]
						special_out = special_special_out[idx] #list Query GCF
						if not args.redundant:
							if special_out!='#':
								asset=set()
								asset.add(special_out.split('|')[1])
								queryDict[query[0]+'#'+str(q)]=asset
							else:
								if exceptionalWP_out[:3]=='WP_':
									accession_from_wp_exceptional=exceptional_wp_results.get(exceptionalWP_out, False)
									if accession_from_wp_exceptional:
										asset=set()
										for elements in accession_from_wp_exceptional:
											asset.add(elements)
										if len(asset)>0:
											if query[0]!=exceptionalWP_out:
												queryDict[exceptionalWP_out+'#'+str(q)+'.'+query[0]]=sortGCFvsGCA(asset)
											else:
												queryDict[exceptionalWP_out+'#'+str(q)]=sortGCFvsGCA(asset)
									else:
										ne+=1
										print(query[0], file= fbad)
						if args.redundant:
							if exceptionalWP_out[1:3]=='P_':
								accession_from_wp_exceptional=exceptional_wp_results.get(exceptionalWP_out, False)
								if accession_from_wp_exceptional:
									asset=set()
									for elements in accession_from_wp_exceptional:
										asset.add(elements)
									if len(asset)>0:
										if query[0]!=exceptionalWP_out:
											queryDict[exceptionalWP_out+'#'+str(q)+'.'+query[0]]=sortGCFvsGCA(asset)
										else:
											queryDict[exceptionalWP_out+'#'+str(q)]=sortGCFvsGCA(asset)
								else:
									ne+=1
									print(query[0], file= fbad)
							elif exceptionalWP_out[2]!='_':
								assembly_from_identical=identicalProtID_redundant(identicalProtID_Out)
								if assembly_from_identical!='#':
									asset=set()
									for elements in assembly_from_identical:
										asset.add(elements)
									if len(asset)>0:
										#queryDict[identicalProtID_Out+'#'+str(q)+'.'+query[0]]=sortGCFvsGCA(asset)
										if query[0]!=exceptionalWP_out:
											queryDict[exceptionalWP_out+'#'+str(q)+'.'+query[0]]=sortGCFvsGCA(asset)
										else:
											queryDict[exceptionalWP_out+'#'+str(q)]=sortGCFvsGCA(asset)
								else:
									ne+=1
									print(query[0], file= fbad)

			else:
				if query[0][:3]=='XP_' and query[0][-2]=='.': #XP Accession
					asset=set()
					accession_from_xp_out=xp_paired_results.get(query[0], False)
					if accession_from_xp_out:
						for bioprojs in accession_from_xp_out:
							if bioprojs in bioDict:
								if bioDict[bioprojs]==query[1]:
									asset.add(query[1])
						if len(asset)>0:
							queryDict[query[0]+'#'+str(q)]=asset
					else:
						ne+=1
						print(query[0], file= fbad)
				elif query[0][:3]!='XP_' and query[0][-2]=='.': #not XP Accession
					asset=set()
					accession_from_wp_out=wp_paired_results.get(query[0], False)
					if accession_from_wp_out:
						for elements in accession_from_wp_out:
							if query[1]==elements:
								asset.add(query[1])
					if len(asset)>0:
						queryDict[query[0]+'#'+str(q)]=asset
				else:
					ne+=1
					print(query[0], file= fbad)
else:
	for query in queryList:
		q+=1
		queryDict[query[0]+'#'+str(q)]=query[1]

nai=0
NqueryDict={} #{'WP_019504790.1#1': ['GCF_000332195.1'], 'WP_028108719.1#2': ['GCF_000422645.1'], 'WP_087820443.1#3': ['GCF_900185565.1']}
if args.localGenomeList:
	with open (args.out_prefix+'_Insufficient_Info_In_DB.txt', 'w') as fNai:
		for query in queryDict:
			assemblyIdlist=[]
			if queryDict[query]!={'NAI'}:
				assemblyId=queryDict[query]
				faa_gz=localDir+queryDict[query]+'.faa.gz'
				if os.path.isfile(faa_gz):
					with gzip.open(localDir+assemblyId+'.faa.gz', 'rb') as faaIn:
						for line in faaIn:
							if line.decode('utf-8')[0]=='>':
								Line=line.decode('utf-8').rstrip()
								if '>'+query.split('#')[0]==Line.split(' ')[0]:
									gff_gz=localDir+assemblyId+'.gff.gz'
									if os.path.isfile(gff_gz):
										with gzip.open(localDir+assemblyId+'.gff.gz', 'rb') as lgffIn: #Download and read gff.gz
											cds_c=0
											name_c=0
											prot_c=0
											cset=set()
											nset=set()
											pset=set()
											for line in lgffIn:
												if line.decode('utf-8')[0]!='#':
													Line=line.decode('utf-8').rstrip().split('\t')
													if Line[2]=='CDS':
														cds_c=1
														cset.add(cds_c)
														if Line[8].split(';')[3][:5]=='Name=': #eliminates pseudo gene as they don't have 'Name='
															name_c=1
															nset.add(name_c)
															if Line[8].split(';')[3].split('=')[1]==query.split('#')[0]:
																assemblyIdlist.append(assemblyId)
																NqueryDict[query]=list(set(assemblyIdlist))
																prot_c=1
																pset.add(prot_c)
															else:
																prot_c=0
																pset.add(prot_c)
														else:
															name_c=0
															nset.add(name_c)
													else:
														cds_c=0
														cset.add(cds_c)
											if lcheck(pset)>0:
												pass
											elif lcheck(pset)==0:
												if lcheck(cset)>0 and lcheck(nset)>0:
													print(query.split('#')[0],' did not match with supplement local GFF File')
													print(query, file=fNai)
													nai+=1
											else:
												print('Use recommended [NCBI refseq] format of GFF file')
												break
									else:
										print("Error: %s file not found" % gff_gz)
				else:
					print("Error: %s file not found" % faa_gz)
			else:
				print(query, file=fNai)
				nai+=1
else:
	with open (args.out_prefix+'_Insufficient_Info_In_DB.txt', 'w') as fNai:
		for query in queryDict:
			if not queryDict[query] or len(queryDict[query])==0 or queryDict[query]=={'NAI'}:
				print(query, file=fNai)
				nai+=1
				continue
			if args.redundant:
				redun=0
				for newRed in (redundantCreate(queryDict[query],args.redundant)):
					redun+=1
					NqueryDict[query+'.'+str(redun)]=list(str(newRed).split())
			else:
				NqueryDict[query]=random.sample(queryDict[query],1)

#print(NqueryDict)

if not args.localGenomeList:
	print('\n> Downloading Genome Assembly Files from NCBI FTP Server \n')

_print_lock = threading.Lock()

speciesNameFromOnlineDict = {}

if not args.localGenomeList:
	# Build flat list of (query, item, 1-based-index) jobs
	jobs = []
	for query, items in NqueryDict.items():
		for item in items:
			jobs.append((query, item))

	# Worker count: respect -c flag, cap at 8 to avoid hammering NCBI
	max_workers = min(int(args.cpu) if args.cpu else 8, len(jobs), 8)

	completed_count = 0
	with ThreadPoolExecutor(max_workers=max_workers) as executor:
		future_to_job = {
			executor.submit(_download_assembly, query, item, idx + 1, len(jobs)): (query, item)
			for idx, (query, item) in enumerate(jobs)
		}
		for future in as_completed(future_to_job):
			item_result, species_label, success = future.result()
			completed_count += 1
			if species_label is not None:
				speciesNameFromOnlineDict[item_result] = species_label

newQ = len(NqueryDict)

for query in NqueryDict:
	for item in NqueryDict[query]:
		if args.localGenomeList:
			faaFile=item+'.faa.gz'
			fastaSeq = gzip.open(localDir+faaFile, "rt")
			for record in SeqIO.parse(fastaSeq, "fasta"):
				if record.id==query.split('#')[0]:
					if args.redundant:
						speciesNameFromOnlineDict[item]=remBadChar(record.description.split('[')[-1][:-1])+'_'+remBadChar(item)
					else:
						speciesNameFromOnlineDict[item]=remBadChar(record.description.split('[')[-1][:-1])
		else:
			if item not in speciesNameFromOnlineDict or speciesNameFromOnlineDict[item]=='Nothing':
				faaFile=item+'.faa.gz'
				if os.path.isfile(localDir+faaFile):
					fastaSeq = gzip.open(localDir+faaFile, "rt")
					for record in SeqIO.parse(fastaSeq, "fasta"):
						if record.id==query.split('#')[0]:
							if args.redundant:
								speciesNameFromOnlineDict[item]=remBadChar(record.description.split('[')[-1][:-1])+'_'+remBadChar(item)
							else:
								speciesNameFromOnlineDict[item]=remBadChar(record.description.split('[')[-1][:-1])
				else:
					speciesNameFromOnlineDict[item]=item+'#not_found'


if args.keep:
	with open(args.out_prefix+'_speciesInfo.txt','w') as asmOut:
		for query in NqueryDict:
			for item in NqueryDict[query]:
				print(item, query.split('#')[0], speciesNameFromOnlineDict[item], sep='\t', file=asmOut)

#print(NqueryDict) #{'WP_019504790.1#1': ['GCF_000332195.1'], 'WP_028108719.1#2': ['GCF_000422645.1'], 'WP_087820443.1#3': ['GCF_900185565.1']}

print('\n'+'>> Input file assessment report: ')
print('\t'+'Discarded protein ids with improper accession : '+str(ne)+'. See "'+args.out_prefix+'_NameError.txt'+'" file for details.')
print('\t'+'Discarded protein ids lacking proper information in RefSeq DB : '+str(nai)+'. See "'+args.out_prefix+'_Insufficient_Info_In_DB.txt'+'" file for details.')
print('\t'+'Remaining queries: '+str(newQ))

FoundDict={} #Accession that found in Refseq
FlankFoundDict={} #Accession that have flanking genes
accFlankDict={} #{'WP_092250023.1#1': {0: 'WP_092250023.1+', 1: 'WP_092250020.1+', 2: 'WP_092250017.1-', -1: 'tRNA*+', -2: 'WP_092250026.1-'}}
positionDict={} #Accession as keys:Start and end position as value
speciesDict={} #SpeciesName stored here
queryStrand={} #Strand Information for each query
LengthDict={} #Length of each query
seqDict={}
desDict={}
acc_CGF_Dict={}
treeFastadict={} #Query as key and sequence in fasta as value
querySeqDict={} #For Tree Command


count=0
for query in NqueryDict:
	count+=1
	if args.verbose:
		print('\n'+'> '+str(count)+' in process out of '+str(newQ)+' ... '+'\n')
		print('Query Name =', query.split('#')[0], '\n')
	for item in NqueryDict[query]:  #{'WP_019504790.1#1': ['GCF_000332195.1'],NqueryDict
		a=0
		LineList=[]
		geneProt={} # 'gene2504': 'WP_092248795.1', 'gene1943': 'tRNA'
		geneChrom={} #'gene1708': 'NZ_MJLP01000030.1'
		#item='GCF_'+items[4:]
		speciesNameFromDB=speciesNameFromOnlineDict[item]
		gff_gz=localDir+item+'.gff.gz'
		if os.path.isfile(gff_gz):
			LineList=[]
			geneProt={} # 'gene2504': 'WP_092248795.1', 'gene1943': 'tRNA'
			geneChrom={} #'gene1708': 'NZ_MJLP01000030.1'
			with gzip.open(localDir+item+'.gff.gz', 'rb') as gffIn: #Download and read gff.gz
				for line in gffIn:
					if line.decode('utf-8')[0]!='#':
						Line=line.decode('utf-8').rstrip().split('\t')
						if Line[2]=='CDS':
							if Line[8].split(';')[3][:5]=='Name=': #eliminates pseudo gene as they don't have 'Name='
								if query.split('_')[0]=='XP':
									if 'GeneID:' in Line[8]:
										geneProt[getGeneId(Line[8])]=Line[8].split(';')[3].split('=')[1]
										geneChrom[getGeneId(Line[8])]=Line[0]
										#print(getGeneId(Line[8]), Line[8].split(';')[3].split('=')[1], Line[0], 'gpc#')
								else:
									#print(Line[8].split(';')[1].split('=')[1], Line[8].split(';')[3].split('=')[1], Line[0], 'gpc#')
									geneProt[Line[8].split(';')[1].split('=')[1]]=Line[8].split(';')[3].split('=')[1]
									geneChrom[Line[8].split(';')[1].split('=')[1]]=Line[0]
									##print(geneProt)>'GeneID:187667': 'NP_493855.2' NP_417570.1 gene-b3099
						if Line[2][-4:]=='gene':
							a+=1
							if query.split('_')[0]=='XP':
								if 'GeneID:' in Line[8]:
									newGene=str(a)+'\t'+getGeneId_gene(Line[8])+'\t'+ Line[3]+'\t'+Line[4]+'\t'+ Line[6]+ '\t'+ Line[0]
									#print(newGene) #1	GeneID:353377	3747	3909	-	NC_003279.8
									LineList.append(newGene.split('\t'))
									for genDes in Line[8].split(';'):
										if 'gene_biotype=' in genDes:
											if getGeneId_gene(Line[8]) not in geneProt:
												geneProt[getGeneId_gene(Line[8])]=genDes.split('=')[1]+'_'+query.split('#')[1]+'.'+str(random.randint(0,int(s)*2-1))+'*'
							else:
								newGene=str(a)+'\t'+Line[8].split(';')[0][3:]+'\t'+ Line[3]+'\t'+Line[4]+'\t'+ Line[6]+ '\t'+ Line[0]
								LineList.append(newGene.split('\t')) #1	   gene3006		10266   10342   -	   NZ_FPCC01000034.1
								for genDes in Line[8].split(';'):
									if 'gene_biotype=' in genDes:
										if Line[8].split(';')[0][3:] not in geneProt:
											geneProt[Line[8].split(';')[0][3:]]=genDes.split('=')[1]+'_'+query.split('#')[1]+'.'+str(random.randint(0,int(s)*2-1))+'*'
				geneList=[] ##List of gene names coding same protein (accession), we are taking one from them
				for genes in geneProt:
					if geneProt[genes]==query.split('#')[0]:
						geneList.append(genes)
				if len(geneList)>0:
					rangeList=[]
					for line in LineList:
						if geneChrom[geneList[0]]==line[5]:
							rangeList.append(int(line[0]))
					for genes in geneProt:
						if genes==geneList[0]:
							if query.split('#')[0]==geneProt[genes]:
								for line in LineList:
									if genes==line[1]:
										FoundDict[query]='Yes'
										if args.tree or args.hmmdb:
											treeFastadict[query]=str(seqFasLocal(item,query))
											querySeqDict[query+'|'+remBadChar(spLocal(item,query.split('#')[0]))]=str(seqLocal(item,query.split('#')[0]))
										if speciesNameFromDB!='Nothing' or speciesNameFromDB!='':
											#print(speciesNameFromDB, 1)
											speciesDict[query]=speciesNameFromDB
										else:
											#print(spLocal(item, query.split('#')[0]), 2)
											speciesDict[query]=spLocal(item, query.split('#')[0])
										lineIdx = LineList.index(line)
										query_num = query.split('#')[1]
										queryStrand[query]= LineList[lineIdx][4]
										positionDict[query]= ("\t".join(map(str,LineList[lineIdx][2:-2])))
										LengthDict[query]= int(LineList[lineIdx][3])-int(LineList[lineIdx][2])+1
										udsDict={}
										dsDict={}
										udsDict[0]= query+'+' # O strand
										lengthCheck=[]
										for x in range(1,int(s)):
											if lineIdx-x>=0 and lineIdx-x<len(LineList):
												if rangeList.count(int(LineList[lineIdx-x][0]))>0:
													acc_CGF_Dict[query]= LineList[lineIdx-x][-1] +'\t'+ item
													seqDict[str(geneProt[LineList[lineIdx-x][1]])]=localNone(seqLocal(item, geneProt[LineList[lineIdx-x][1]]))
													desDict[geneProt[LineList[lineIdx-x][1]]]=desLocal(item, geneProt[LineList[lineIdx-x][1]])
													positionDict[geneProt[LineList[lineIdx-x][1]]+'#'+query_num]= ("\t".join(map(str,LineList[lineIdx-x][2:-2])))
													lengthCheck.append(seqFasLenLocal(item,geneProt[LineList[lineIdx-x][1]]+'#'+query_num))
													LengthDict[geneProt[LineList[lineIdx-x][1]]+'#'+query_num]= int(LineList[lineIdx-x][3])-int(LineList[lineIdx-x][2])+1
													udsDict[int(ups(LineList[lineIdx][4])[0]+str(x))]= geneProt[LineList[lineIdx-x][1]]+'#'+query_num+\
														normalize_strand(LineList[lineIdx][4],LineList[lineIdx-x][4])
										for y in range(1,int(s)):
											if lineIdx+y<len(LineList):
												if rangeList.count(int(LineList[lineIdx+y][0]))>0:
													acc_CGF_Dict[query]= LineList[lineIdx+y][-1] +'\t'+ item
													seqDict[str(geneProt[LineList[lineIdx+y][1]])]=localNone(seqLocal(item,geneProt[LineList[lineIdx+y][1]]))
													desDict[geneProt[LineList[lineIdx+y][1]]]=desLocal(item,geneProt[LineList[lineIdx+y][1]])
													positionDict[geneProt[LineList[lineIdx+y][1]]+'#'+query_num]= ("\t".join(map(str,LineList[lineIdx+y][2:-2])))
													lengthCheck.append(seqFasLenLocal(item,geneProt[LineList[lineIdx+y][1]]+'#'+query_num))
													LengthDict[geneProt[LineList[lineIdx+y][1]]+'#'+query_num]= int(LineList[lineIdx+y][3])-int(LineList[lineIdx+y][2])+1
													dsDict[int(downs(LineList[lineIdx][4])[0]+str(y))]= geneProt[LineList[lineIdx+y][1]]+'#'+query_num+\
														normalize_strand(LineList[lineIdx][4],LineList[lineIdx+y][4])
										if lenChecker(lengthCheck)=='keep':
											udsDict.update(dsDict)
											accFlankDict[query]=udsDict
											if query in accFlankDict:
												if len(accFlankDict[query])>0:
													FlankFoundDict[query]='Yes'
													if args.verbose:
														print('\t', query.split('#')[0], 'Report: Flanking Genes Found', '\n')
												else:
													FlankFoundDict[query]='No'
													if args.verbose:
														print('\t', query.split('#')[0], 'Report: Flanking Genes Not Found', '\n')
											else:
												FlankFoundDict[query]='No'
												if args.verbose:
													print('\t', query.split('#')[0], 'Report: Flanking Genes Not Found', '\n')
										else:
											FlankFoundDict[query]='Yes'
											if args.verbose:
												print('\t', query.split('#')[0], 'Report: Flanking Gene is longer than 5000 amino acid, thus query is discarded', '\n')
								else:
									if query not in FoundDict:
										FlankFoundDict[query]='No'
										if args.verbose:
											print('\t', query.split('#')[0], 'Report: Flanking Genes Not Found', '\n')

				else:
					FlankFoundDict[query]='No'
					FoundDict[query]='No: ProteinID was not found in Genome Assembly'
					if args.verbose:
						print('\t', query.split('#')[0], 'Report: Flanking Genes Not Found', '\n')
		else:
			FlankFoundDict[query]='No'
			FoundDict[query]='No: ProteinID was not found in Genome Assembly'
			if args.verbose:
				print('\t', query.split('#')[0], 'Report: Flanking Genes Not Found', '\n')

if not args.localGenomeList:
	if args.keep:
		pass
	else:
		subprocess.Popen("rm GC*_*.gz", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

allFlankGeneList=[]
for keys in accFlankDict:
	for item in accFlankDict[keys]:
		allFlankGeneList.append(accFlankDict[keys][item].split('#')[0])


flankF=0
with open (args.out_prefix+'_flankgene_Report.log', 'w') as errOut:
	serial=0
	print('#Serial','Query','Assembly_Found', 'FlankingGene_Found', sep='\t', file=errOut)
	for queries in NqueryDict:
		serial+=1
		if queries in FoundDict:
			if queries in FlankFoundDict:
				if FlankFoundDict[queries]=='Yes':
					flankF+=1
				print(str(serial), queries.split('#')[0], FoundDict[queries], FlankFoundDict[queries], sep='\t', file=errOut)
			else:
				print(str(serial), queries.split('#')[0], FoundDict[queries], 'No', sep='\t', file=errOut)
		else:
			print(str(serial), queries.split('#')[0], 'No', 'No', sep='\t', file=errOut)

reportDict={}
for query in queryList:
	queryNumList=[]
	for queryNum in NqueryDict:
		if query[0] in queryNum:
			queryNumList.append(queryNum)
	if len(queryNumList)>0:
		reportDict[query[0]]=queryNumList
	else:
		reportDict[query[0]]=str('no').split()

qcount=0
discardedGene=0
with open (args.out_prefix+'_QueryStatus.txt', 'w') as sumOut:
	print('#Serial', 'Status', sep='\t', file=sumOut)
	for query in queryList:
		#print('query',query)
		qcount+=1
		for item in reportDict[query[0]]:
			#print('item',item)
			if item in FlankFoundDict:
				if item in accFlankDict:
					if FlankFoundDict[item]=='Yes':
						print(str(qcount), reporter(query[0], item.split('#')[0], similarityID(query[0], item.split('#')[0]), ''.join(NqueryDict[item]), 'Yes'), sep='\t', file=sumOut)
					else:
						print(str(qcount), reporter(query[0], item.split('#')[0], similarityID(query[0], item.split('#')[0]), ''.join(NqueryDict[item]), 'No'), sep='\t', file=sumOut)
				else:
					if FlankFoundDict[item]=='Yes':
						discardedGene+=1
						print(str(qcount), query[0]+' is a valid NCBI protein accession but Discarded : Flanking gene with a length more than 5000 amino acid detected.', sep='\t',file=sumOut)
					else:
						print(str(qcount), reporter(query[0], item.split('#')[0], similarityID(query[0], item.split('#')[0]), ''.join(NqueryDict[item]), 'No'), sep='\t', file=sumOut)
			else:
				print(str(qcount), reporter(query[0], 'No', 'No', 'No', 'No'), sep='\t', file=sumOut)

print('\n'+'>> Flanking Genes found : '+str(flankF)+' out of remaining '+str(serial)+'. See "'+args.out_prefix+'_flankgene_Report.log'+'" file for details.'+'\n'+'\n')

if int(flankF)==0:
	print('>> No Flanking Genes found, please update your accession list.')
	sys.exit()
elif len(accFlankDict)==0:
	print('>> For every query, flanking gene(s) having a length more than 5000 amino acid detected. Please update your accession list. \n')
	sys.exit()
elif int(flankF)==discardedGene:
	print('>> For every query, flanking gene(s) having a length more than 5000 amino acid detected. Please update your accession list. \n')
	sys.exit()
else:
	pass


if args.tree: #Generate fasta file for making phylogenetic Tree
	with open(args.out_prefix+'_tree.fasta', 'w') as treeOut:
		for queries in NqueryDict:
			if queries in accFlankDict:
				if queries in treeFastadict:
					print(treeFastadict[queries], file=treeOut)


if len(seqDict)!=len(desDict):
	if len(seqDict)>len(desDict):
		for seqids in sorted(seqDict):
			if seqids not in desDict:
				desDict[seqids]=des_check(str(seq_from_wp(seqids).split('\t')[0]))
	else:
		for seqids in sorted(desDict):
			if seqids not in seqDict:
				seqDict[seqids]=str(seq_from_wp(seqids).split('\t')[1])
else:
	if args.verbose:
		print ('Description collected for Flanking Genes!')

with open(args.out_prefix+'_all.fasta', 'w') as all_fasta:
	for queries in NqueryDict:
		if queries in accFlankDict:
			if queries in treeFastadict:
				print(treeFastadict[queries], file=all_fasta)
	for seqids in sorted(seqDict):
		if seqDict[seqids]!='--':
			print('>'+desDict[seqids]+'\n'+seqDict[seqids], file=all_fasta)


b=0
with open (args.out_prefix+'_flankgene.fasta'+'_cluster_out', 'w') as fastaNew:
	for seqids in sorted(seqDict):
		if seqDict[seqids]!='--':
			b+=1
			print('>'+seqids+'|'+desDict[seqids]+'\n'+seqDict[seqids], file=fastaNew)

if args.verbose:
	print ('Total Flanking genes found = '+ str(b))

print('\n>> Now running Jackhmmer and clustering flanking genes\n')

infilename=args.out_prefix+'_flankgene.fasta'+'_cluster_out'

directory = args.out_prefix+'_flankgene.fasta'+'_cluster_out_individuals'
if not os.path.exists(directory):
	os.makedirs(directory)

infile=open(args.out_prefix+'_flankgene.fasta'+'_cluster_out',"r").read()
al=infilename+"_"+iters+"_"+evthresh+"_jackhits.tsv"
outacclists=open(al,"w")

percentileJack=0
i=1
for seqids in sorted(seqDict): #Running Jackhmmer for finding homologs
	if seqDict[seqids]!='--':
		percentileJack+=1
		if args.verbose:
			if percentileJack % 5 == 0:
				print('\t'+'>>> '+str(round(int(percentileJack)*100/b))+'%'+' Completed...'+'('+str(percentileJack)+'/'+str(b)+')')
			if percentileJack % b == 0:
				print('\t'+'>>> Completed ' +'\n')
		i_f=directory+"/"+str(i)+".txt"
		indivfile=open(i_f,"w")
		indivfile.write(">"+seqids+'\n'+seqDict[seqids])
		indivfile.close()
		if args.cpu:
			command="jackhmmer --cpu %s -N %s --incE %s --incdomE %s --tblout %s/tblout%s.txt %s  %s>%s/out%s.txt" %(core, iters, evthresh, evthresh, directory, str(i), i_f, infilename, directory, str(i))
		else:
			command="jackhmmer -N %s --incE %s --incdomE %s --tblout %s/tblout%s.txt %s  %s>%s/out%s.txt" %(iters, evthresh, evthresh, directory, str(i), i_f, infilename, directory, str(i))
		os.system(command)
		tbl=open(directory+"/tblout"+str(i)+".txt").read()
		part=tbl.split("----------\n")[1].split("\n#")[0]
		lines=part.splitlines()
		acclist=[]
		for line in lines:
			lineList=line.split()
			if len(lineList)>17:
				inc=line.split()[17]
				acc=line.split('|')[0]
				if inc=="1":
					acclist.append(acc)
		outacclists.write(str(i)+"\t"+str(acclist)+"\n")
		i=i+1

outacclists.close()

raw=open(al).read().strip()

d={}
index=0
for line in raw.split("\n"):
	if line.split("\t")[1]!='[]':
		index+=1
		actxt=line.split("\t")[1].replace(",","").replace("[","").replace("]","").replace("'","")
		actlist=actxt.split(" ")
		d[index]=(actlist)

i=1
while i<len(d)+1:	#use i and j to iterate through the combinations
	list1=d[i]
	j=i+1
	while j<len(d)+1:
		list2=d[j]
		#print "i", i, list1, " vs ",
		#print "j", j, list2
		if set(list1) & (set(list2)): # if there is an intersection
			#print "yes there is", list1, list2
			union=list(set(list2) | set(list1))
			d[j]=union
			d[i]=[] #...and empty the redundant list
		j=j+1
	i=i+1

#d : {1: ['WP_000153877.1'], 2: [], 3: ['WP_000291520.1'], 4: ['WP_000342211.1'], 5: [],... 30: ['WP_055032751.1', 'WP_001246052.1']}

trueAccessionCount={}
for keys in d:
	numbers=[]
	for item in d[keys]:
		numbers.append(allFlankGeneList.count(item))
	trueAccessionCount[(';'.join(map(str,d[keys])))]=sum(numbers)

#trueAccessionCount: 'WP_001229260.1;WP_001229255.1': 4
odtrue=OrderedDict(sorted(trueAccessionCount.items(), key= lambda item:item[1],reverse=True))

familyNumber=0
with open(infilename+"_"+iters+"_"+evthresh+"_clusters.tsv","w") as clusOut:
	for k, v in odtrue.items():
		#print (k.split(';'), odtrue[k], v)
		if len(k.split(';'))>0 and v>0:
			familyNumber+=1
			print(str(familyNumber),str(odtrue[k]),k, sep='\t', file=clusOut)

outfile_des=open(infilename+"_"+iters+"_"+evthresh+"_outdesc.txt","w")
inf=open(infilename+"_"+iters+"_"+evthresh+"_clusters.tsv","r")

acclists=inf.read().splitlines()
for line in acclists:
	acclist=line.split("\t")[2].split(";")
	familyAssignedValue=line.split("\t")[0]
	if int(line.split("\t")[1])>1:
		for acc in acclist:
			outfile_des.write(familyAssignedValue+'('+str(allFlankGeneList.count(acc))+')'+"\t"+acc+"\t"+desDict[acc]+"\n")
		outfile_des.write ("\n\n")

outfile_des.close()

#domain search using hmmscan with default output file + domain table format - under construction
if args.hmmdb:
	print('\n>> Now running hmmscan and searching domains\n')
	if args.cpu:
		hmmscan_cmd = "hmmscan -E 1e-10 --cpu %s -o %s_dom.txt --domtblout %s_dom_out.txt %s %s_all.fasta"%(round(core/3), args.out_prefix, args.out_prefix, args.hmmdb, args.out_prefix)
	else:
		hmmscan_cmd = "hmmscan -E 1e-10 -o %s_dom.txt --domtblout %s_dom_out.txt %s %s_all.fasta"%(args.out_prefix, args.out_prefix, args.hmmdb, args.out_prefix)
	subprocess.run(hmmscan_cmd, shell=True) #runs hmmscan

	#domain search ends
	#domain dictionary for outdesc output
	dom_dict={}
	dom_key=''
	with open(args.out_prefix+"_dom_out.txt",'r') as dom_dict_infile:
		for line in dom_dict_infile:
			if any(e in line for e in seqDict.keys()):
				dom_dict_line = line.split()
				if dom_dict_line[3]!=dom_key:
					dom_key = dom_dict_line[3]
					dom_key_cut=dom_key.split('|')
					dom_value = "\t"+"\t"+dom_dict_line[0]+"\t"+dom_dict_line[5]+"\t"+dom_dict_line[6]+"\t"+dom_dict_line[17]+"\t"+dom_dict_line[18]+"\n"
					dom_dict[dom_key_cut[0]]=dom_value


	with open(infilename+"_"+iters+"_"+evthresh+"_outdesc.txt", 'r') as out_desc_in, open(infilename+"_"+iters+"_"+evthresh+"_outdesc_dom.txt", 'w') as out_desc_out:
		for line in out_desc_in:
			if any(e in line for e in dom_dict.keys()):
				split_line = line.split('\t')
				out_desc_out.write(line)
				out_desc_out.write(''.join(dom_dict[split_line[1]]))
			else:
				out_desc_out.write(line)

familyDict={} # Accession:Assigned family Number from Jackhammer
with open(args.out_prefix+'_flankgene.fasta_cluster_out_'+iters+'_'+evthresh+'_clusters.tsv', 'r') as clusterIn:
	for line in clusterIn:
		if line[0]!='#':
			line=line.rstrip().split('\t')
			if int(line[1])>1:
				for item in (line[2].split(';')):
					familyDict[item]=int(line[0])
			else:
				familyDict[line[2]]=0

familynum=[]
for acc in familyDict:
	familynum.append(familyDict[acc])

center=int(max(familynum))+1
noProt=int(max(familynum))+2
noProtP=int(max(familynum))+3
noColor=int(max(familynum))+4
for ids in LengthDict:
	if ids.split('#')[0][-1]=='*':
		if ids.split('#')[0][:2].lower()=='ps':
			familyDict[ids.split('#')[0]]=noProtP
		else:
			familyDict[ids.split('#')[0]]=noProt
	if ids.split('#')[0] not in familyDict:
		if ids in NqueryDict:
			familyDict[ids.split('#')[0]]=center
		else:
			familyDict[ids.split('#')[0]]=noColor

color={}
color[noColor]='#ffffff'
color[center]='#000000'
color[noProt]='#f2f2f2'
color[noProtP]='#f2f2f3'
colorDict={}  #Assigned family Number from Jackhammer : colorcode
for families in set(familynum):
	if families == 0:
		colorDict[families]=str('#ffffff')
	else:
		if random_color()!='#ffffff' or random_color()!='#000000' or random_color()!='#f2f2f2' or random_color()!='#f2f2f3' :
			colorDict[families]=random_color()

colorDict.update(color)

maxs=(int(s)-1) # required to calculate border size of postscript output
mins=maxs-(maxs*2) # required to calculate border size of postscript output

if not args.tree_order:
	write_operon_tsv(args.out_prefix+'_operon.tsv', accFlankDict)
	draw_operon_pdf(args.out_prefix+'_operon.tsv', args.out_prefix+'_operon.pdf')

if args.tree:###Tree Command with ETE###
	tree_file= args.out_prefix+'_tree.fasta'
	if args.cpu:
		tree_command="ete3 build -a %s -o %s --nochecks --clearall -w mafft_default-trimal01-none-fasttree_full --rename-dup-seqnames --cpu %s" %(tree_file, tree_file[:-6], core)
	else:
		tree_command="ete3 build -a %s -o %s --nochecks --clearall -w mafft_default-trimal01-none-fasttree_full --rename-dup-seqnames" %(tree_file, tree_file[:-6])
	os.system(tree_command)
	from ete3 import Tree, SeqMotifFace, TreeStyle, add_face_to_node

	def normalize_strandView(item):  #Strand view change
		if item=='+':
			return '>'
		else:
			return '<'

	def familyView(item):  #Strand view change
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
			return str(item)

	seqMult=((maxs)*2)+1
	seq = ("XXXXXXXXXXXXX--"*seqMult)
	startDict={}
	udList=[]
	for ud in range (mins, maxs+1, 1):
		udList.append(ud)
	sList=[]
	for sa in range(1, 15*seqMult, 15):
		sList.append(sa)
	for ln in range(len(udList)):
		startDict[udList[ln]]=sList[ln]

	nwTree=''
	motifDict={}
	motifDict_2={}
	if os.path.isfile(args.out_prefix+'_tree/mafft_default-trimal01-none-fasttree_full/'+args.out_prefix+'_tree.fasta.final_tree.nw') == True:
		with open(args.out_prefix+'_tree/mafft_default-trimal01-none-fasttree_full/'+args.out_prefix+'_tree.fasta.final_tree.nw', 'r') as treeIn:
			for line in treeIn:
				nwTree=line
				for items in line.replace('(','').replace(')', '').replace(';', '').replace(',','\t').split('\t'):
					item=items.split('|')[0]
					simple_motifs=[]
					simple_motifs_2=[]
					for keys in sorted(startDict):
						if keys in accFlankDict[item]:
							simple_motifs_s = [startDict[keys], startDict[keys]+13, normalize_strandView(accFlankDict[item][keys][-1]), None, size, outliner(colorDict[familyDict[accFlankDict[item][keys][:-1].split('#')[0]]]), 'rgradient:'+colorDict[familyDict[accFlankDict[item][keys][:-1].split('#')[0]]], "arial|"+fsize+"|black|"+familyView(familyDict[accFlankDict[item][keys][:-1].split('#')[0]])]
							simple_motifs.append(simple_motifs_s)
							simple_motifs_2_s = [startDict[keys], startDict[keys]+13, normalize_strandView(accFlankDict[item][keys][-1]), None, size, outliner(colorDict[familyDict[accFlankDict[item][keys][:-1].split('#')[0]]]),colorDict[familyDict[accFlankDict[item][keys][:-1].split('#')[0]]], "arial|"+fsize+"|black|"]
							simple_motifs_2.append(simple_motifs_2_s)
						else:
							simple_motifs_s = [startDict[keys], startDict[keys]+13, '[]', None, size, '#eeeeee', 'rgradient:'+'#ffffff', "arial|"+fsize+"|black|"]
							simple_motifs.append(simple_motifs_s)
							simple_motifs_2_s = [startDict[keys], startDict[keys]+13, '[]', None, size, '#eeeeee', '#ffffff', "arial|"+fsize+"|black|"]
							simple_motifs_2.append(simple_motifs_2_s)
					motifDict[items[:items.index(':')]]=simple_motifs
					motifDict_2[items[:items.index(':')]]=simple_motifs_2
	else:
		print('> ETE3 failed to create tree due to lack of valid protein accesions, at least 2 required !')
		sys.exit()

	def get_example_tree():
		# Create a random tree and add to each leaf a random set of motifs
		# from the original set
		t= Tree(nwTree)
		for item in nwTree.replace('(','').replace(')', '').replace(';', '').replace(',','\t').split('\t'):
			seqFace = SeqMotifFace(seq, motifs=motifDict[item[:item.index(':')]], seq_format="-", gap_format="blank")
			(t & item[:item.index(':')]).add_face(seqFace, 0, "aligned")
		t.ladderize()
		return t

	def get_example_tree_2():
		# Create a random tree and add to each leaf a random set of motifs
		# from the original set
		t= Tree(nwTree)
		for item in nwTree.replace('(','').replace(')', '').replace(';', '').replace(',','\t').split('\t'):
			seqFace2 = SeqMotifFace(seq, motifs=motifDict_2[item[:item.index(':')]], seq_format="-", gap_format="blank")
			(t & item[:item.index(':')]).add_face(seqFace2, 0, "aligned")
		t.ladderize()
		return t

	if __name__ == '__main__':
		t = get_example_tree()
		ts = TreeStyle()
		ts.tree_width = 300
		ts.show_branch_support = True
		if args.tree_order:
			t.write(outfile=args.out_prefix+'_ladderTree.nw')
			t.render(args.out_prefix+"_flankgenes_1.svg",tree_style=ts)
		else:
			t.render(args.out_prefix+"_flankgenes_1.svg",tree_style=ts)


	if __name__ == '__main__':
		t = get_example_tree_2()
		ts = TreeStyle()
		ts.tree_width = 300
		ts.show_branch_support = True
		if args.tree_order:
			t.write(outfile=args.out_prefix+'_ladderTree.nw')
			t.render(args.out_prefix+"_flankgenes_2.svg", tree_style=ts)
		else:
			t.render(args.out_prefix+"_flankgenes_2.svg", tree_style=ts)

if args.tree and args.tree_order:  # Queries in postscript file will be presented as tree order
	treeOrderList=[]
	with open(args.out_prefix+'_ladderTree.nw', 'r') as laddertreeIn:
		for line in laddertreeIn:
			for items in line.replace('(','').replace(')', '').replace(';', '').replace(',','\t').split('\t'):
				item=items.split('|')[0]
				treeOrderList.append(item)
	write_operon_tsv(args.out_prefix+'_TreeOrder_operon.tsv', treeOrderList)
	draw_operon_pdf(args.out_prefix+'_TreeOrder_operon.tsv', args.out_prefix+'_TreeOrder_output.pdf')

print('\n'+'<<< Done >>>')
print('\nIf you use FlaGs2 in your work, please remember to cite these papers!'+'\n\n- Saha CK, Pires RS, Brolin H, Delannoy M, Atkinson GC. 2020. FlaGs and webFlaGs: discovering novel biology through the analysis of gene neighbourhood conservation. Bioinformatics.'+\
'\nhttps://doi.org/10.1093/bioinformatics/btaa788'+\
'\n\n- Jimmy S, Saha CK, Kurata T, Stavropoulos C, Oliveira SRA, Koh A, Cepauskas A, Takada H, Rejman D, Tenson T, Strahl H, Garcia-Pino A, Hauryliuk V, Atkinson GC (2020).\nA widespread toxin-antitoxin system exploiting growth control via alarmone signaling. Proc. Natl. Acad. Sci. U.S.A.'+\
'\nhttps://doi.org/10.1073/pnas.1916617117\n')
endtime = time.perf_counter()
print(f"Executed in: {endtime - starttime:.6f} seconds")
sys.exit()
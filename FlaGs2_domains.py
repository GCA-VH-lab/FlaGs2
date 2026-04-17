# This script runs hmmscan and generates the FlaGs graphical output as an .html file, the content of which can be downloaded as .svg. The output file also offers some interactivity.


import plotly.graph_objs as go
from plotly.subplots import make_subplots
import csv
import colorsys
import random
import argparse
import subprocess
from FlaGs2 import postscriptSize, random_color, _hls2hex, outliner, parse_arguments

args = parse_arguments()

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
newQ=0
#domain search using hmmscan with default output file + domain table format
def hmmscaner(args):
    print('\n>> Now running hmmscan and searching domains\n')
    if args.cpu:
        hmmscan_cmd = "hmmscan -E 1e-10 --cpu %s -o %s_dom.txt --domtblout %s_dom_out.txt %s %s_all.fasta"%(round(core/3), args.out_prefix, args.out_prefix, args.hmmdb, args.out_prefix)
    else:
        hmmscan_cmd = "hmmscan -E 1e-10 -o %s_dom.txt --domtblout %s_dom_out.txt %s %s_all.fasta"%(args.out_prefix, args.out_prefix, args.hmmdb, args.out_prefix)
    subprocess.run(hmmscan_cmd, shell=True)
    dom_dict=[]
    with open(args.out_prefix+"_dom_out.txt",'r') as dom_dict_infile, open(args.out_prefix+"_domains.tsv",'w') as dom_out:
        for line in dom_dict_infile:
            if not line.startswith("#"):
                dom_dict_line = line.split()
                dict_line = []
                dom_key_cut=dom_dict_line[3].split('|')
                dict_line.append(dom_key_cut[0])
                dict_line.append(dom_dict_line[6])
                dict_line.append(dom_dict_line[0])
                dict_line.append(dom_dict_line[19])
                dict_line.append(dom_dict_line[20])
                dom_dict.append(dict_line)
                dom_out.write('\t'.join(dict_line))
                dom_out.write('\n')


hmmscaner(args)
color={}
colorDict={}
colorDict_domains={}

# 1. Merging the two traces into one plot
fig = make_subplots(shared_yaxes = True, shared_xaxes = True)


# 2. Setting the axes. Since the output is a single graph, 
#    the y-axis labels need to be set a specified
with open(args.out_prefix+'_TreeOrder_operon.tsv', newline='') as f:
    rows = list(csv.reader(f, delimiter='\t'))

unique_col0 = len(set(row[0] for row in rows if row))

PX_PER_ROW   = 40         # pixels allocated per operon row
y_axes = unique_col0 * PX_PER_ROW

# x-width: derive from the actual genomic coordinate span across all genes,
# then add space for the left tick labels and the right-side legend.
x_coords = [int(r[5]) for r in rows if len(r) > 6 and r[5].lstrip('-').isdigit()] + [int(r[6]) for r in rows if len(r) > 6 and r[6].lstrip('-').isdigit()]
coord_span = (max(x_coords) - min(x_coords)) if x_coords else 5000
SCALE     = 0.15  # genomic coords → pixels (tune if genes appear too squished)
x_axes = int(coord_span * SCALE)


# 3. All lists
arrowList = []
domainList = []
xList_gene = []
yList_gene = []
xList_domain = []
yList_domain = []
y_tick_marks = []
labels = []
genes = []


# 4. Data file input
main_file = open('example_TreeOrder_operon.tsv','r').read()
eg1 = main_file.split("\n\n\n\n")
y_level_m = 0
for m in eg1:
    if m != '':
        row1 = 0
        entries1 = m.splitlines()
        ndoms=len(entries1)
        y_level_m = y_level_m-10-round(postscriptSize(newQ))
        for entry in entries1:
            items1 = entry.split("\t")
            x_gene_start = int(items1[5])
            x_gene_end = int(items1[6])
            dx_gene_length = int(items1[1])
            gene_direction = items1[3]
            gen1_name = int(items1[4])
            accesssion = str(items1[0])
            id1 = str(items1[9])                     

            # 4a. When genes are to small the arrow shape is distorted because the coordinates are too close to each other.
            #     This makes these genes longer to keep the shape of the arrow. 
            if dx_gene_length < 100:
                x_gene_start = int(items1[5])-50
                x_gene_end = int(items1[6])+50
            else:
                x_gene_start = int(items1[5])
                x_gene_end = int(items1[6])


            #Colours (imported from FlaGs script)
            center=int(gen1_name)+1
            noProt=int(gen1_name)+2
            noProtP=int(gen1_name)+3
            noColor=int(gen1_name)+4
            
            color[noColor]='#ffffff'
            color[center]='#000000'
            color[noProt]='##f2f2f2'
            color[noProtP]='#f2f2f3'


            if gen1_name == 0:
                colorDict[gen1_name]=str('#ffffff')
            elif x_gene_start == 1:                    
                colorDict[gen1_name]=str('#000000')
            elif 'pseudogene_' in id1:
                colorDict[gen1_name]=str('#f2f2f2') 
            elif 'tRNA_' in id1:
                colorDict[gen1_name]=str('#f2f2f3')
            else:
                if gen1_name not in colorDict:
                    colorDict[gen1_name] = random_color()


            # 4b. Editing the label for each gene in the legend
            protein = (str(id1) + ' ' + '(' + ('Start: {}\tEnd: {}'.format(x_gene_start, x_gene_end)) + ')')
            hover_text = 'ID: ' + id1 +'<br>Start: ' + str(x_gene_start) + '<br>End: ' + str(x_gene_end)


            # 4c. Drawing the genes as polygons/arrows
            if gene_direction == '-':
                xList_gene = [x_gene_start+100, x_gene_start, x_gene_start+100, x_gene_end, x_gene_end, x_gene_start+100]
                yList_gene = [y_level_m-2, y_level_m, y_level_m+2, y_level_m+2, y_level_m-2, y_level_m-2]
                arrowList.append(fig.add_trace(go.Scatter(x = xList_gene, y = yList_gene, fill="toself", fillcolor='rgba(0,0,0,0)', opacity = 1, line=dict(color = outliner(colorDict[gen1_name])), mode = 'lines+text', name = protein, showlegend = False)))
            else:
                xList_gene = [x_gene_start, x_gene_start, x_gene_end-100, x_gene_end, x_gene_end-100, x_gene_start]
                yList_gene = [y_level_m-2, y_level_m+2, y_level_m+2, y_level_m, y_level_m-2, y_level_m-2] 
                arrowList.append(fig.add_trace(go.Scatter(x = xList_gene, y = yList_gene, fill="toself", fillcolor='rgba(0,0,0,0)', opacity = 1, line=dict(color = outliner(colorDict[gen1_name])), mode = 'lines+text', name = protein, showlegend = False)))
            



  
            # 6. Domains file input
            domain_file = open(args.out_prefix+"_domains.tsv", 'r').read()
            eg2 = domain_file.split("\n")
            id1 = str(items1[9])
            for d in eg2:
                if d != '':
                    row2 = 0
                    entries2=d.splitlines()
                    ndoms=len(entries2)
                    y_level_d = y_level_m 
                    for entry2 in entries2:
                        if entry2 == '':
                            continue # go to end of loop
                        items2 = entry2.split('\t')
                        x_domain_start = int(items2[3]) + x_gene_start
                        x_domain_end = int(items2[4]) + x_gene_start
                        dx_domain_size = int(x_domain_end)-int(x_domain_start)
                        id2 = str(items2[0])
                        domain_name = str(items2[2])

                        # 6a. Editing the label for each domain in the legend
                        domain = ('     ' + id2 + ' ' + '(' + ('Start: {}\tEnd: {}'.format(x_domain_start, x_domain_end)) + ')')
                        domain_group = (domain_name)

                        # 6b. If a gene has additional information about domains (i.e. same id is found in second file), then these will also be drawn inside the arrow.
                        if id1.startswith(id2):
                            if domain_name not in colorDict_domains:
                                colorDict_domains[domain_name] = random_color()
                            if x_domain_end != x_gene_end and x_domain_start != x_gene_start:
                                xList_domain = [x_domain_start, x_domain_start, x_domain_end, x_domain_end, x_domain_start]
                                yList_domain = [y_level_m-2, y_level_m+2, y_level_m+2, y_level_m-2, y_level_m-2]
                                domainList.append(fig.add_trace(go.Scatter(x=xList_domain, y=yList_domain, fill="toself", hoverinfo = 'none', fillcolor=colorDict_domains[domain_name], line=dict(color=colorDict_domains[domain_name]), opacity = 0.5, mode='lines', name = domain, legendgroup=domain_name,legendgrouptitle_text=domain_group)))        
                            elif x_domain_end == x_gene_end and gene_direction == '-':
                                xList_domain = [x_domain_start, x_domain_start, x_domain_end, x_domain_end, x_domain_start]
                                yList_domain = [y_level_m-2, y_level_m+2, y_level_m+2, y_level_m-2, y_level_m-2]
                                domainList.append(fig.add_trace(go.Scatter(x=xList_domain, y=yList_domain, fill="toself", hoverinfo = 'none', fillcolor=colorDict_domains[domain_name], line=dict(color=colorDict_domains[domain_name]), opacity = 0.5, mode='lines', name = domain, legendgroup=domain_name,legendgrouptitle_text=domain_group)))
                            elif x_domain_start == x_gene_start and gene_direction == '-':
                                xList_domain = [x_domain_start+100, x_domain_start, x_domain_start+100, x_domain_end, x_domain_end, x_domain_start+100]
                                yList_domain = [y_level_m-2, y_level_m, y_level_m+2, y_level_m+2, y_level_m-2, y_level_m-2]
                                domainList.append(fig.add_trace(go.Scatter(x=xList_domain, y=yList_domain, fill="toself", hoverinfo = 'none', fillcolor=colorDict_domains[domain_name], line=dict(color=colorDict_domains[domain_name]), opacity = 0.5, mode='lines', name = domain, legendgroup=domain_name,legendgrouptitle_text=domain_group)))
                            elif x_domain_start == x_gene_start and gene_direction == '+':
                                xList_domain = [x_domain_start, x_domain_start, x_domain_end, x_domain_end, x_domain_start]
                                yList_domain = [y_level_m-2, y_level_m+2, y_level_m+2, y_level_m-2, y_level_m-2]
                                domainList.append(fig.add_trace(go.Scatter(x=xList_domain, y=yList_domain, fill="toself", hoverinfo = 'none', fillcolor=colorDict_domains[domain_name], line=dict(color=colorDict_domains[domain_name]), opacity = 0.5, mode='lines', name = domain, legendgroup=domain_name,legendgrouptitle_text=domain_group)))                                       
                            elif x_domain_end == x_gene_end and gene_direction == '+':
                                xList_domain = [x_domain_start, x_domain_start, x_domain_end-100, x_domain_end, x_domain_end-100, x_domain_start]
                                yList_domain = [y_level_m-2, y_level_m+2, y_level_m+2, y_level_m, y_level_m-2, y_level_m-2]
                                domainList.append(fig.add_trace(go.Scatter(x=xList_domain, y=yList_domain, fill="toself", hoverinfo = 'none', fillcolor=colorDict_domains[domain_name], line=dict(color=colorDict_domains[domain_name]), opacity = 0.5, mode='lines', name = domain, legendgroup=domain_name,legendgrouptitle_text=domain_group)))                                       
                            else:
                                pass


            # 7. Setting the y labels i.e. the organism name and accession nr etc.
            y_tick_marks += [y_level_m]
            labels += [items1[0]]
            genes += [items1[9]]

            row1 = row1+1


# 8. Changing the download format of the .html as .svg instead of the defaul .png
config = {'toImageButtonOptions': {'format': 'svg','filename': 'FlaGs','scale': 1}}


# 9. Graph layout
fig.update_xaxes(visible = False)
fig.update_yaxes(visible = True, showgrid = False, showline = False, autorange = True, automargin = True, showticklabels = True, tickvals = y_tick_marks, ticktext = labels, ticklen = 20, tickmode = 'array', titlefont = dict(family = 'Open Sans', size = 8))
fig.update_layout(autosize=False, width=x_axes, height=y_axes, margin=dict(b=10, t=10, pad=10), paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)', showlegend = True)

fig.show(config=config)

fig.write_html(args.out_prefix+"domain_output.html")
# This script runs hmmscan and generates the FlaGs graphical output as an .html file, the content of which can be downloaded as .svg. The output file also offers some interactivity.


import plotly.graph_objs as go
from plotly.subplots import make_subplots
import csv
import colorsys
import random
import argparse
import subprocess
from FlaGs2 import postscriptSize, _hls2hex, outliner, parse_arguments, random_color

args = parse_arguments()

newQ = 0

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

def contrast_color(color_number):
    r=3
    s=7
    c = (color_number%r/r) + color_number//r/(s+color_number//(r*s)*2)
    if c < 0: c = 1 - c
    if c > 1: c = c - 1
    d = 0.5
    e = 0.5
    return _hls2hex(c, d, e)

# 1. Merging the two traces into one plot
fig = make_subplots(shared_yaxes = True, shared_xaxes = True)


# 2. Setting the axes. Since the output is a single graph, 
#    the y-axis labels need to be set a specified
with open(args.out_prefix+'_TreeOrder_operon.tsv', newline='') as f:
    rows = list(csv.reader(f, delimiter='\t'))

unique_col0 = len(set(row[0] for row in rows if row))

PX_PER_ROW  = 40
y_axes = unique_col0 * PX_PER_ROW
x_coords = [int(r[5]) for r in rows if len(r) > 6 and r[5].lstrip('-').isdigit()] + [int(r[6]) for r in rows if len(r) > 6 and r[6].lstrip('-').isdigit()]
coord_span = (max(x_coords) - min(x_coords)) if x_coords else 5000
SCALE     = 0.12  # genomic coords to pixels
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
unique_domain_names = []
eg1 = []
eg2 = []

color={}
colorDict_domains={}
colorDict={}

# Data file input
main_file = open('example_TreeOrder_operon.tsv','r').read()
eg1 = main_file.split("\n\n\n\n")
# Domain file input
with open (args.out_prefix+"_domains.tsv", 'r') as domain_file:
    for line in domain_file:
        eg2.append(line.split('\t'))
eg2.sort(key = lambda x: int(x[3]))
y_level_m = 0
for m in eg1:
    if m != '':
        row1 = 0
        entries1 = m.splitlines()
        y_level_m = y_level_m-10-round(postscriptSize(newQ))
        for entry in entries1:
            prev_labels = []
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

            # 4c. Drawing the genes as polygons/arrows
            if gene_direction == '-':
                xList_gene = [x_gene_start+100, x_gene_start, x_gene_start+100, x_gene_end, x_gene_end, x_gene_start+100]
                yList_gene = [y_level_m-2, y_level_m, y_level_m+2, y_level_m+2, y_level_m-2, y_level_m-2]
            else:
                xList_gene = [x_gene_start, x_gene_start, x_gene_end-100, x_gene_end, x_gene_end-100, x_gene_start]
                yList_gene = [y_level_m-2, y_level_m+2, y_level_m+2, y_level_m, y_level_m-2, y_level_m-2] 
            arrowList.append(fig.add_trace(go.Scatter(
                x = xList_gene, 
                y = yList_gene, 
                fill="toself", 
                fillcolor='rgba(0,0,0,0)', 
                opacity = 1, 
                line=dict(color = outliner(colorDict[gen1_name]), width = 1), 
                mode = 'lines+text', 
                showlegend = False, 
                hoverinfo = 'none')))
            
            for items2 in eg2:
                        id2 = str(items2[0])
                        if id1.startswith(id2):
                            if gene_direction == '-':
                                x_domain_start = x_gene_end - int(items2[3])
                                x_domain_end = x_gene_end - int(items2[4])
                                direction_multiplier = -1
                            else:
                                x_domain_start = int(items2[3]) + x_gene_start
                                x_domain_end = int(items2[4]) + x_gene_start
                                direction_multiplier = 1
                            domain_name = str(items2[2])
                            domain = ('     ' + id2 + ' ' + '(' + ('Start: {}\tEnd: {}'.format(x_domain_start, x_domain_end)) + ')')
                            isfirst = False # check if it is first domain of this kind for legend grouping
                            if not domain_name in unique_domain_names:
                                isfirst = True
                                unique_domain_names.append(domain_name)
                                colorDict_domains[domain_name] = contrast_color(len(colorDict_domains))
                            domain_number = str(unique_domain_names.index(domain_name))
                            domain_group = domain_number + '. ' + domain_name
                            xList_domain = [x_domain_start, x_domain_start, x_domain_end, x_domain_end, x_domain_start]
                            yList_domain = [y_level_m-2, y_level_m-1, y_level_m+2, y_level_m-2, y_level_m-2]
                            domainList.append(fig.add_trace(go.Scatter(
                                x=xList_domain, 
                                y=yList_domain, 
                                fill="toself", 
                                fillcolor=colorDict_domains[domain_name], 
                                opacity = 1, 
                                line=dict(color='rgba(0,0,0,0)'), 
                                mode='lines', 
                                name = domain_group, 
                                legendgroup=domain_name, 
                                showlegend=isfirst)))
                            # Label: domain numbers above the containing protein
                            if any(abs(prev_label - (x_domain_start+x_domain_end)/2)<175 for prev_label in prev_labels):
                                label_start=prev_labels[-1]+175*direction_multiplier
                            else:
                                label_start = (x_domain_start+x_domain_end)/2
                            domainList.append(fig.add_trace(go.Scatter(
                                x=[label_start],
                                y=[y_level_m + 2],
                                mode='text',
                                text=[domain_number],
                                textposition='top center',
                                textfont=dict(family='Open Sans', size=12, color=colorDict_domains[domain_name]),
                                hoverinfo='none',
                                showlegend=False,
                                legendgroup=domain_name,)))
                            prev_labels.append(label_start)
                            if gene_direction == '-':
                                prev_labels.sort(reverse=True)
                            else:
                                prev_labels.sort()

            # 7. Setting the y labels i.e. the organism name and accession nr etc.
            y_tick_marks += [y_level_m]
            raw = items1[0]; parts = raw.split('_', 3); labels += ['_'.join(parts[:3])]

            row1 = row1+1


# 8. Changing the download format of the .html as .svg instead of the defaul .png
config = {'toImageButtonOptions': {'format': 'svg','filename': 'FlaGs','scale': 1}}


# 9. Graph layout
fig.update_xaxes(visible = False)
fig.update_yaxes(visible = True, showgrid = False, showline = False, autorange = True, automargin = True, showticklabels = True, tickvals = y_tick_marks, ticktext = labels, ticklen = 20, tickmode = 'array', titlefont = dict(family = 'Open Sans', size = 8))
fig.update_layout(autosize=False, width=x_axes, height=y_axes, margin=dict(b=10, t=10, pad=10), paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)', legend=dict(orientation="h",   yanchor="top",  y=0, xanchor="left", x=0), showlegend = True)

fig.show(config=config)

fig.write_html(args.out_prefix+"domain_output.html")
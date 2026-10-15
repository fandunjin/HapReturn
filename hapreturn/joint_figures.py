"""Joint AGP/PAF figures: paired coordinate evidence and paired-cell source matrix.

Tracks use saved interval evidence, not smoothed or expanded window labels.
Matrix halves use the same whole-target-length denominator. Descriptive data
from one assembly; uncertainty intervals and hypothesis tests do not apply.
"""
import csv
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch, Rectangle
from matplotlib.colors import LinearSegmentedColormap

STYLE = {'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
         'font.size': 7, 'axes.titlesize': 8, 'axes.labelsize': 7, 'xtick.labelsize': 6,
         'ytick.labelsize': 6, 'svg.fonttype': 'none', 'pdf.fonttype': 42,
         'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': .5,
         'legend.frameon': False, 'savefig.facecolor': 'white'}


def read_rows(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def export_figure(fig, out_dir, stem, axes):
    fig.canvas.draw()
    # Record measured geometry for portable report generation.
    boxes = [a.get_window_extent().transformed(fig.dpi_scale_trans.inverted()) for a in axes]
    if len(boxes) == 2:
        for key in ('y0', 'y1', 'width', 'height'):
            if abs(getattr(boxes[0], key)-getattr(boxes[1], key))*72 > 1.5:
                raise ValueError('Comparable figure panels are misaligned')
    (out_dir/(stem+'.layout.json')).write_text(json.dumps({'units':'inches', 'axes':[list(b.bounds) for b in boxes]}, indent=2))
    audit_dir = os.environ.get('NATURE_FIGURE_AUDIT_DIR')
    if audit_dir:
        sys.path.insert(0, audit_dir)
        from audit_panel_alignment import require_matplotlib_panel_alignment
        require_matplotlib_panel_alignment(fig, json_out=out_dir/(stem+'.alignment.json'),
                                           tolerance_pt=1.5, strict=True,
                                           exclude_axes=[a for a in fig.axes if a not in axes])
    fig.savefig(out_dir/(stem+'.svg'))
    fig.savefig(out_dir/(stem+'.pdf'))
    fig.savefig(out_dir/(stem+'.png'), dpi=600)
    plt.close(fig)


def make_joint_figures(out_dir, labels, rows, segments):
    out_dir = Path(out_dir)
    palette = {labels[0]:'#397ca7', labels[1]:'#d99a50', 'unresolved':'#c4cbd2', 'both':'#8e83a6'}
    names = [r['seq_id'] for r in rows]
    by_name = {n: [s for s in segments if s['seq_id']==n] for n in names}
    groups = [[r for r in rows if r['assignment']==l] for l in labels]
    unresolved = [r for r in rows if r['assignment'] not in labels]
    # Keep unresolved records visible in the second column with an explicit label.
    groups[1] += unresolved
    count = max(len(g) for g in groups)
    maximum = max(r['length'] for r in rows)/1e6
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, 2, figsize=(183/25.4, max(150, count*9+40)/25.4), squeeze=True)
        fig.subplots_adjust(left=.14, right=.98, bottom=.11, top=.87, wspace=.44)
        for panel, (ax, group, label) in enumerate(zip(axes, groups, labels)):
            ax.set_title(label + (' + unresolved' if panel==1 and unresolved else ''), loc='left', pad=10)
            ax.annotate(chr(97+panel), (0,1), xycoords='axes fraction', xytext=(-38,10),
                        textcoords='offset points', fontsize=9, fontweight='bold')
            for i, row in enumerate(group):
                name = row['seq_id']
                for lane, offset, height in [('A',-.24,.16), ('M',0,.16), ('H',.24,.23)]:
                    for s in by_name[name]:
                        if lane=='A':
                            source=s['agp_source']
                        elif lane=='M':
                            support=[x for x in s['alignment_sources'].split(',') if x]
                            source=support[0] if len(support)==1 else 'both' if len(support)>1 else 'unresolved'
                        else:
                            source=s['source']
                        ax.add_patch(Rectangle(((s['start']-1)/1e6,i+offset-height/2),s['length_bp']/1e6,height,
                                               facecolor=palette.get(source,palette['unresolved']), edgecolor='none'))
                        if lane=='H' and s['conflict']:
                            ax.add_patch(Rectangle(((s['start']-1)/1e6,i+offset-height/2),s['length_bp']/1e6,height,
                                                   fill=False,edgecolor='#b44d5e',linewidth=.6))
                if all(row[l+'_support_bp']>0 for l in labels):
                    ax.plot(row['length']/1e6+.6,i,'o',ms=2,color='#263b50',clip_on=False)
            ax.set_xlim(0,maximum*1.04)
            ax.set_ylim(count-.5,-.6)
            ax.set_yticks(range(len(group)), [r['seq_id'] for r in group])
            for tick, row in zip(ax.get_yticklabels(),group):
                tick.set_gid('select-chr-'+str(names.index(row['seq_id'])))
            ax.set_xlabel('Chromosome position (Mb)')
            ax.tick_params(axis='y',length=0,pad=5)
            ax.spines['left'].set_visible(False)
            ax.grid(axis='x',color='#e9edf1',linewidth=.35)
            ax.set_axisbelow(True)
        fig.legend(handles=[Patch(facecolor=palette[l],label=l) for l in labels]+
                   [Patch(facecolor=palette['unresolved'],label='Unknown / no accepted alignment'),
                    Patch(facecolor=palette['both'],label='Both alignment sources')],
                   loc='upper center',bbox_to_anchor=(.52,.98),ncol=2,fontsize=6)
        fig.text(.14,.045,'Within each chromosome, top to bottom: AGP (A), minimap2 (M), integrated evidence (H).',fontsize=6)
        fig.text(.14,.022,'Shared coordinate scale; white = gap; dot = mixed integrated origin; red outline = conflict.',fontsize=6)
        export_figure(fig,out_dir,'joint_evidence_tracks',list(axes))
    matrix_path=out_dir/'alignment_matrix.tsv'
    component_path=out_dir/'curation/agp_chromosome_components.tsv'
    if not matrix_path.exists() or not component_path.exists():
        return
    length_by_name={r['seq_id']:r['length'] for r in rows}
    origin=defaultdict(float)
    keys=set()
    for p in read_rows(component_path):
        if p['seq_id'] in names:
            key=(p['source_hap'],p['source_contig'])
            keys.add(key)
            origin[(key,p['seq_id'])]+=int(p['length_bp'])/length_by_name[p['seq_id']]
    aligned={}
    for p in read_rows(matrix_path):
        key=(p['source'],p['original_contig'])
        if any(float(p.get(n,0))>0 for n in names):
            keys.add(key)
            for n in names:
                aligned[(key,n)]=float(p.get(n,0))
    order=[r['seq_id'] for group in groups for r in group]
    def row_order(key):
        best=max(order, key=lambda n:max(origin.get((key,n),0), aligned.get((key,n),0)))
        return labels.index(key[0]) if key[0] in labels else 2, order.index(best), key[1]
    keys=sorted(keys,key=row_order)
    # Interleave AGP and PAF fractions in each contig-by-chromosome cell.
    values=np.zeros((len(keys),len(order)*2))
    data=[]
    for i,key in enumerate(keys):
        for j,name in enumerate(order):
            a,m=origin.get((key,name),0),aligned.get((key,name),0)
            if not 0<=a<=1+1e-8 or not 0<=m<=1+1e-8:
                raise ValueError('Evidence matrix fraction outside [0,1]')
            values[i,2*j:2*j+2]=a,m
            data.append({'source':key[0],'original_contig':key[1],'seq_id':name,'agp_fraction':a,'alignment_fraction':m})
    from .agp import write_tsv
    write_tsv(out_dir/'joint_evidence_matrix.tsv',data,list(data[0]) if data else
              ['source','original_contig','seq_id','agp_fraction','alignment_fraction'])
    with plt.rc_context(STYLE):
        fig,ax=plt.subplots(figsize=(183/25.4,max(170,len(keys)*3+44)/25.4))
        fig.subplots_adjust(left=.26,right=.94,top=.90,bottom=.20)
        cmap=LinearSegmentedColormap.from_list('support',['#ffffff','#cadce9','#76a6c4','#245477'])
        image=ax.imshow(values,aspect='auto',vmin=0,vmax=1,cmap=cmap,interpolation='nearest')
        ax.set_xticks([2*i+.5 for i in range(len(order))],order,rotation=90,ha='right',rotation_mode='anchor')
        ax.set_yticks(range(len(keys)),[k[1].removeprefix('lcl|') for k in keys])
        for tick,key in zip(ax.get_yticklabels(),keys):
            tick.set_color(palette.get(key[0],palette['unresolved']))
        ax.tick_params(length=0,pad=3)
        for x in range(1,len(order)):
            ax.axvline(2*x-.5,color='#cbd4dd',lw=.2)
        boundary=len(groups[0])*2-.5
        if 0<boundary<values.shape[1]-.5:
            ax.axvline(boundary,color='#334e68',lw=.8)
        ax.set_title('Joint contig-to-chromosome evidence',loc='left',pad=12)
        ax.set_ylabel('Original hifiasm contig',labelpad=8)
        cbax=fig.add_axes([.50,.085,.30,.012])
        cb=fig.colorbar(image,cax=cbax,orientation='horizontal',ticks=[0,.5,1])
        cb.set_label('Fraction of target chromosome length',fontsize=6)
        cb.ax.tick_params(labelsize=6,length=2)
        fig.text(.26,.962,'Each cell: AGP contribution (left) | minimap2 coverage (right)',fontsize=7)
        fig.text(.26,.025,'Columns follow the recovered partition. Contig label colours indicate original source.',fontsize=6)
        export_figure(fig,out_dir,'joint_evidence_matrix',[ax])

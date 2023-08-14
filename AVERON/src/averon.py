#!/usr/bin/env python
# coding: utf-8
import glob
import json
import ipywidgets as widgets
from ipywidgets import IntProgress, Layout, VBox, HBox
from IPython.display import display

import lifelines
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import median_survival_times
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import requests
import scikit_posthocs as sp
import scipy.stats as ss
from scipy.stats import binom_test
import seaborn as sns
from sklearn import preprocessing
import statsmodels.api as sm
import statsmodels.sandbox.stats.multicomp as mltc
from statsmodels.formula.api import ols
from statsmodels.stats.outliers_influence import OLSInfluence
import tabulate
import tkinter;
import warnings

###Ligands from IUPHAR####
def get_ligands_iuphar(gene):
    import pygtop
    N_approved = 0
    N_not_approved = 0
    out = []
    a = pygtop.get_targets_by({'geneSymbol':gene})
    approved_ligands = []
    all_ligands = []
    if len(a)>0:
        target = a[0]
        ligids = []
        todrop = []
        ligands = target.ligands(species='Human')
        for i in range(0,len(ligands)):
            try:
                ligid = ligands[i].ligand_id()
                if ligid not in ligids:
                    ligids.append(ligid)
                else:
                    todrop.append(i)
            except:
                pass
        for i in sorted(todrop, reverse=True):
            del ligands[i]

        if len(ligands)>0:
            for lig in ligands:
                try:
                    lig_name = lig.name()
                    lig_id = lig.ligand_id()
                    lig_type = lig.ligand_type()
                except:
                    print(lig)
                    lig_name = -1
                if lig_name != -1:
                    if lig.approved():
                        lig_status = "approved"
                    else:
                        lig_status = "not_approved"
                    out.append([lig_name,lig_id,lig_type,lig_status])
                
    return(out)

def show_ligands(out1,progressbar,genes):
    
    import tabulate
    from IPython.display import display, Markdown
    
    ligands_out = []
    noligands = []
    for gene in genes:
        
        ligout = get_ligands_iuphar(gene)
        if len(ligout)==0:
            noligands.append(gene)
            nl1 = [gene,0,0,'','']
        else:
            num_approved = 0
            num_total = 0
            approved = []
            all_ligands = []
            with out1:
                link = widgets.HTML(
                    value='<a href="https://www.genenames.org/tools/search/#!/?query='+gene+'" style="color:#ff0000;" target="_blank">'+gene+'</a>',
                )
                display(link)
            lig_tbl = []
            N = 0
            nl = []
            for line in ligout:
                
                N+=1
                lig_name = line[0]
                all_ligands.append(lig_name)
                lig_id = str(line[1])
                lig_status = line[3]
                num_total+=1
                if lig_status == "approved":
                    num_approved+=1
                    approved.append(lig_name)
                    s1 = '<a href="https://www.guidetopharmacology.org/GRAC/LigandDisplayForward?ligandId='+lig_id+'" style="color:#ff8200;">'+lig_name+'</a>'
                else:
                    s1 = '<a href="https://www.guidetopharmacology.org/GRAC/LigandDisplayForward?ligandId='+lig_id+'">'+lig_name+'</a>'
                nl.append(s1)
                if N >= 5:
                    lig_tbl.append(nl)
                    N = 0
                    nl = []
            nl1 = [gene,num_approved,num_total,(':').join(approved),(':').join(all_ligands)]
            lig_tbl.append(nl)
            lig_tbl = tabulate.tabulate(lig_tbl, tablefmt='unsafehtml')
        
            with out1:
                display(lig_tbl)
        
        progressbar.value += 1
        ligands_out.append(nl1)
        
        
    if len(noligands)>0:
        with out1:
            display(Markdown("<span style='color:red'>No ligands: </span>"))
        nolig_tbl = []
        N = 0
        nl = []   
        for gene in noligands:
            N+=1
            s1 = '<a href="https://www.genenames.org/tools/search/#!/?query='+gene+'" target="_blank">'+gene+'</a>'
            nl.append(s1)
            if N >= 10:
                nolig_tbl.append(nl)
                N = 0
                nl = []
        nolig_tbl.append(nl)
        nolig_tbl = tabulate.tabulate(nolig_tbl, tablefmt='unsafehtml')
        with out1:
            display(nolig_tbl)
    ligands_df = pd.DataFrame(ligands_out,columns = ['gene','num_approved_drugs','num_total_ligands','approved_drugs','all_ligands'])
    ligands_df.set_index('gene',inplace=True)
    return(ligands_df)
###END Ligands from IUPHAR

def get_txt_description_of_regulations(geneA,mutA,partner,pdf,net_dict):
    descr_base = geneA + " " + (";").join(mutA) + "/"+partner + " interaction can regulate "
    out = []
    for i in pdf.index:
        pathway = pdf.loc[i]['SET']
        description = descr_base + pathway
        net_dict["nodes"].append({"data":{"id":pathway,"classes":"pathway"}})
        included_genes = pdf.loc[i]['sig_genes_included']
        included_genes = included_genes.split(";")
        included_genes = [x.strip() for x in included_genes]
        included_genes = (",").join(included_genes)
        description = description + " by upregulating "+ included_genes
        print(description)
        out.append(description)
    return ("\n").join(out)
def create_interactive_network(geneA, partner,enrichment_df,fqval,colors,sig_genes,net_folder):
    '''
    geneA - partner
    partner - genes
    enrichment_df - dataframe generated by pathway_enrichment
    colors: list of colors for mutant, partner, signature genes, and pathways
    net_folder: a path to the output folder
    '''
    import ipycytoscape
    import json
    import ipywidgets
    import json
    net_dict = {"nodes":[],"edges":[]}
    net_dict["nodes"].append({"data":{"id":geneA,'label': geneA, "classes":"mutant"}})
    net_dict["nodes"].append({"data":{"id":partner,'label': partner, "classes":"partner"}})
    net_dict["edges"].append({"data":{"source": geneA, "target": partner}})
    
     
    
    for gene in sig_genes:
        net_dict["nodes"].append({"data":{"id":gene,'label': gene,"classes":"gene"}})
        net_dict["edges"].append({"data":{"source": partner, "target": gene}})

    

    pdf = enrichment_df.loc[enrichment_df['fqval']<fqval]
    pdf.reset_index(inplace=True)
    pdf.drop(['index'],inplace=True,axis=1)

    for i in pdf.index:
        pathway = pdf.loc[i]['SET']
        net_dict["nodes"].append({"data":{"id":pathway,'label': pathway,"classes":"pathway"}})
        included_genes = pdf.loc[i]['sig_genes_included']
        included_genes = included_genes.split(";")
        included_genes = [x.strip() for x in included_genes]
        
        for gene in included_genes:
            if gene not in [geneA, partner]:
                net_dict["edges"].append({"data":{"source": gene, "target": pathway}})

    my_style = [
        {'selector': 'node','style': {
            'font-family': 'arial',
            'font-size': '10px',
            'label': 'data(label)',}},
    
        {'selector': 'node[classes="mutant"]','style': {
            'background-color': colors[0]}},
        {'selector': 'node[classes="partner"]','style': {
            'background-color': colors[1]}},
        {'selector': 'node[classes="gene"]','style': {
            'background-color': colors[2]}},
        {'selector': 'node[classes="pathway"]','style': {
            'background-color': colors[3],
            'text-valign': 'bottom',
            'text-halign': 'center'
        
        }}
    ]
    with open(net_folder+'net_dict.txt', 'w') as convert_file:
        convert_file.write(json.dumps(net_dict))
    ipycytoscape_obj = ipycytoscape.CytoscapeWidget()
    ipycytoscape_obj.graph.add_graph_from_json(net_folder+'net_dict.txt')
    ipycytoscape_obj.set_tooltip_source('id')
    
    ipycytoscape_obj.set_layout(name='breadthfirst',roots=[geneA],animate = True)
    ipycytoscape_obj.set_style(my_style)
    return (ipycytoscape_obj, pdf, net_dict)


def create_interactive_network2(net_df,type_df,colors=['green','red','blue','orange']):
    
    import ipycytoscape
    import json
    import ipywidgets
    import json
    net_dict = {"nodes":[],"edges":[]}
   

    for i in net_df.index:
        node1 = net_df.loc[i]['Partner']
        node2 = net_df.loc[i]['Pathway']
        type1 = type_df.loc[type_df['node']==node1]['type'].values[0]
        type2 = type_df.loc[type_df['node']==node2]['type'].values[0]
        
        net_dict["nodes"].append({"data":{"id":node1,"classes":type1}})
        net_dict["nodes"].append({"data":{"id":node2,"classes":type2}})
        net_dict["edges"].append({"data":{"source": node1, "target": node2}})
    
    my_style = [
        {'selector': 'node[classes="Driver"]','style': {
            'background-color': colors[0]}},
        {'selector': 'node[classes="Partner"]','style': {
            'background-color': colors[1]}},
        {'selector': 'node[classes="pathway"]','style': {
            'background-color': colors[3]}}
    ]
    with open('net_dict1.txt', 'w') as convert_file:
        convert_file.write(json.dumps(net_dict))
    ipycytoscape_obj = ipycytoscape.CytoscapeWidget()
    ipycytoscape_obj.graph.add_graph_from_json('net_dict1.txt')
    ipycytoscape_obj.set_tooltip_source('id')
    driver = type_df.loc[type_df['type']=='Driver']['node'].values[0]
    ipycytoscape_obj.set_layout(name='breadthfirst',roots=[driver],animate = True)
    ipycytoscape_obj.set_style(my_style)
    return (ipycytoscape_obj)



def display_cynetwork2(net_df,type_df, colors=['#FF0000','#00FF00','#0000FF','#FFA500'],layout='force-directed'):
    
    import py4cytoscape as p4c
    from py4cytoscape import gen_node_color_map
    from py4cytoscape import gen_node_size_map
    from py4cytoscape import scheme_c_number_continuous

    
    nodes_arr = []
    edges_arr = []
    for i in net_df.index:
        node1 = net_df.loc[i]['Partner']
        node2 = net_df.loc[i]['Pathway']
        type1 = type_df.loc[type_df['node']==node1]['type'].values[0]
        type2 = type_df.loc[type_df['node']==node2]['type'].values[0]
        
        nodes_arr.append([node1,type1])
        nodes_arr.append([node2,type2])
        edges_arr.append([node1,node2,'interacts'])
    
    
    nodes_df = pd.DataFrame(nodes_arr,columns = ['id','type'])
    edges_df = pd.DataFrame(edges_arr,columns = ['source','target','interaction'])

    p4c.create_network_from_data_frames(nodes=nodes_df, edges=edges_df, title="network1")
    #colors = ['#FF0000','#00FF00','#0000FF','#FFA500']
    p4c.set_node_color_mapping('type',['Driver','Partner','Pathway'],colors,'d')
    p4c.set_node_shape_default('ELLIPSE')
    p4c.set_node_size_default(20)
    p4c.analyze_network()
    p4c.set_node_size_mapping(**gen_node_size_map('Degree', scheme_c_number_continuous(25, 60)))
    p4c.layout_network(layout_name = layout)
    return p4c


def display_cynetwork(geneA,partner,cancer,sig_genes,enrichment_df,fqval,
                      colors=['#FF0000','#00FF00','#0000FF','#FFA500'],layout='force-directed'):
    '''
    Creat a network neoPPI-signature genes-pathways
    geneA - partner
    partner - genes
    genes - pathways
    colors: list of colors for 'mutant','partner','sig_gene',and 'pathway'
    layout: 'attribute-grid', 'spherical', 'degree-circle', 'box',
             'attributes-layout', 'kamada-kawai', 'force-directed',
             'grid3D', 'cose', 'flatten', 'hierarchical', 'center3d',
             'attribute-circle', 'fruchterman-rheingold-3D', 
             'stacked-node-layout', 'circular', 'grid',
             'fruchterman-rheingold', 'isom', 'force-directed-cl'
    more info: https://github.com/cytoscape/cytoscape-automation/blob/master/for-scripters/Python/advanced-cancer-networks-and-data.ipynb
    '''
    import py4cytoscape as p4c
    from py4cytoscape import gen_node_color_map
    from py4cytoscape import gen_node_size_map
    from py4cytoscape import scheme_c_number_continuous

    
    nodes_arr = []
    edges_arr = []
    nodes_arr.append([geneA,'mutant'])
    nodes_arr.append([partner,'partner'])
    edges_arr.append([geneA,partner,'interacts'])

    for gene in sig_genes:
        if gene not in [geneA,partner]:
            nodes_arr.append([gene,'sig_gene'])
            edges_arr.append([partner,gene,'regulate'])

    pdf = enrichment_df.loc[enrichment_df['fqval']<fqval]
    pdf.reset_index(inplace=True)
    pdf.drop(['index'],inplace=True,axis=1)

    for i in pdf.index:
        pathway = pdf.loc[i]['SET']
        nodes_arr.append([pathway,'pathway'])
        included_genes = pdf.loc[i]['sig_genes_included']
        included_genes = included_genes.split(";")
        included_genes = [x.strip() for x in included_genes]
        for gene in included_genes:
            edges_arr.append([gene,pathway,'pathway'])

    nodes_df = pd.DataFrame(nodes_arr,columns = ['id','type'])
    edges_df = pd.DataFrame(edges_arr,columns = ['source','target','interaction'])

    p4c.create_network_from_data_frames(nodes=nodes_df, edges=edges_df, title=(" ").join([cancer,geneA,partner]))
    #colors = ['#FF0000','#00FF00','#0000FF','#FFA500']
    p4c.set_node_color_mapping('type',['mutant','partner','sig_gene','pathway'],colors,'d')
    p4c.set_node_shape_default('ELLIPSE')
    p4c.set_node_size_default(20)
    p4c.analyze_network()
    p4c.set_node_size_mapping(**gen_node_size_map('Degree', scheme_c_number_continuous(10, 30)))
    p4c.layout_network(layout_name = layout)
    return p4c
    
    

def func1():
    def showit(slider1,out):
        with out:
            print(slider1)
    
    def show_value(slider1):
        showit(slider1,out)

    out = widgets.Output()
    slider1 = widgets.FloatSlider(min=0,max=1,step=0.01,value=0.5)
    widget_1 = widgets.interactive(show_value,slider1 = slider1)    
    vBox_1 = VBox([slider1, out]) 
    widget_1.update()         
    return vBox_1
def get_signature_genes2(correlation_df,partner,geneA,CORR_BP_MUT,PVAL_BP_MUT,PVAL,QVAL):

        sig_genes = correlation_df[partner].loc[(correlation_df[partner]['CORR_BP_MUT']>=CORR_BP_MUT) &
                                (correlation_df[partner]['PVAL_BP_MUT']<=PVAL_BP_MUT) &
                                (correlation_df[partner]['CORR_BP_MUT']>correlation_df[partner]['CORR_BP_WT']) &
                                (correlation_df[partner]['PVAL']<=PVAL) &
                                (correlation_df[partner]['QVAL']<=QVAL)]
        sig_genes = sig_genes.index
    
    
        sig_tbl = []
        N = 0
        nl = []
        for s in sig_genes:
            if N > 9:
                sig_tbl.append(nl)
                N = 0
                nl = []
            s1 = '<a href="https://www.genenames.org/tools/search/#!/?query='+s+'">'+s+'</a>'
            nl.append(s1)
            N+=1
        sig_tbl.append(nl)
        sig_tbl = tabulate.tabulate(sig_tbl, tablefmt='unsafehtml')

        return(sig_genes,sig_tbl)
    
def display_signature_genes(correlation_df,partner,geneA):
    from ipywidgets import Layout
    def get_signature_genes(out,correlation_df,partner,geneA,CORR_BP_MUT,PVAL_BP_MUT,PVAL,QVAL):
        from IPython.display import display, Markdown
        global sig_genes
        
        if CORR_BP_MUT >= 0:
            sig_genes = correlation_df[partner].loc[(correlation_df[partner]['CORR_BP_MUT']>=CORR_BP_MUT) &
                                (correlation_df[partner]['PVAL_BP_MUT']<=PVAL_BP_MUT) &
                                (correlation_df[partner]['CORR_BP_MUT']>correlation_df[partner]['CORR_BP_WT']) &
                                (correlation_df[partner]['PVAL']<=PVAL) &
                                (correlation_df[partner]['QVAL']<=QVAL)]
        else:
            sig_genes = correlation_df[partner].loc[(correlation_df[partner]['CORR_BP_MUT']<=CORR_BP_MUT) &
                                (correlation_df[partner]['PVAL_BP_MUT']<=PVAL_BP_MUT) &
                                (correlation_df[partner]['CORR_BP_MUT']<correlation_df[partner]['CORR_BP_WT']) &
                                (correlation_df[partner]['PVAL']<=PVAL) &
                                (correlation_df[partner]['QVAL']<=QVAL)]
        sig_genes = sig_genes.index
    
    
        sig_tbl = []
        N = 0
        nl = []
        for s in sig_genes:
            if N > 9:
                sig_tbl.append(nl)
                N = 0
                nl = []
            s1 = '<a href="https://www.genenames.org/tools/search/#!/?query='+s+'">'+s+'</a>'
            nl.append(s1)
            N+=1
        sig_tbl.append(nl)
        sig_tbl = tabulate.tabulate(sig_tbl, tablefmt='unsafehtml')
    
        out.clear_output(True)
        header = "\n"+geneA+" MUT/"+partner+" "+str(len(sig_genes))+" signature genes"
        with out:
            display(Markdown("<span style='color:red'>"+header+"</span>"))
            display(sig_tbl)
        return(sig_genes,sig_tbl)
    
    def show_genes(CORR_BP_MUT,PVAL_BP_MUT,PVAL,QVAL):
        global sig_genes
        print("AAA")
        #QVAL=1
        sig_genes, table = get_signature_genes(out,correlation_df,partner,geneA,CORR_BP_MUT,PVAL_BP_MUT,PVAL,QVAL)
    
    CORR_BP_MUT = 0.33
    PVAL_BP_MUT = 0.05
    PVAL = 0.05
    QVAL = 0.25
    
    out = widgets.Output(layout=Layout(height='300px'))

    sliders = widgets.interactive(show_genes,
                    CORR_BP_MUT=widgets.FloatSlider(min=-1,max=1,step=0.01,value=CORR_BP_MUT,
                                        style= {'description_width': 'initial'},
                                        description="CORR_BP_MUT"
                                        ),
                    PVAL_BP_MUT=widgets.FloatSlider(min=0,max=1,step=0.001,value=PVAL_BP_MUT,readout_format = '.3f',
                                        style= {'description_width': 'initial'},
                                        description="PVAL_BP_MUT"),
                    PVAL=widgets.FloatSlider(min=0,max=1,step=0.001,value=PVAL,readout_format = '.3f',
                                 style= {'description_width': 'initial'},
                                        description="PVAL_MUT_vs_WT"),
                    QVAL=widgets.FloatSlider(min=0,max=1,step=0.001,value=QVAL,readout_format = '.3f',
                                 style= {'description_width': 'initial'},
                                        description="QVAL_MUT_vs_WT"))
    
    controls_1 = VBox(sliders.children[:-1])
    vbox = VBox([controls_1,out])
    sliders.update()
    return(vbox,sig_genes)
    
    


def get_outliers_correct(x,y):
    model = sm.OLS(x,y).fit()#formula=f, data=arr).fit()
    threshold = 4/len(x)
    influence = model.get_influence()
    cooks = influence.cooks_distance
    (distance, p_value) = influence.cooks_distance
    infl_data = [i for i in range(0,len(distance)) if distance[i] > threshold]
    return(infl_data)
    
def get_corr_genes(scores,exp_arr,exp_genes):
    out = []
    for i in range(0,len(exp_arr)):
        outliers = get_outliers_correct(scores,exp_arr[i]) #Indices of the outlier points
        new_scores = np.delete(scores,outliers)
        new_exp_vals = np.delete(exp_arr[i],outliers)
        r,pval = ss.pearsonr(new_scores,new_exp_vals)
        out.append([exp_genes[i],r,pval,len(new_scores)])
    return out

def ppi_score_mrna_corr(gene,geneA,partner,all_mut_score_df,all_wt_score_df,df_mut_exp_samples,df_wt_exp_samples,color_mut,color_wt):
    scores_mut = all_mut_score_df[partner].values.tolist()
    gene_exp = df_mut_exp_samples.loc[gene].values.tolist()
    ###Temporary for grant
    todrop = [x for x in range(0,len(gene_exp)) if gene_exp[x]==0]
    scores_mut = [scores_mut[i] for i in range(0,len(scores_mut)) if i not in todrop]
    gene_exp = [gene_exp[i] for i in range(0,len(gene_exp)) if i not in todrop]

    ###Remove it later
    r = 0
    pval = 1
    try:
        if len(set(scores_mut))>3 and len(set(gene_exp))>3:
            #outliers = get_otliers(scores_mut,gene_exp).index.values   
            #outliers = get_outliers1(scores_mut,gene_exp)#.index.values   
            outliers = get_outliers_correct(scores_mut,gene_exp)
            tmp_df = pd.DataFrame()
            tmp_df['scores']=scores_mut
            tmp_df['exp']=gene_exp

            if len(outliers)>0:
                #tmp_df.drop(tmp_df.index[outliers],axis=0,inplace=True)
                tmp_df.drop(outliers,axis=0,inplace=True)
            if len(set(tmp_df['scores'].values.tolist()))>2 and len(set(tmp_df['exp'].values.tolist()))>2:
                r,pval = ss.pearsonr(tmp_df['scores'].values.tolist(),tmp_df['exp'].values.tolist())
            N = len(tmp_df)
    except ValueError:
        r = 0
        pval = 1
        N = 0

    
    #fig, ax = plt.subplots()
    fig = plt.figure(figsize=(10,5))
    grid = plt.GridSpec(12, 2, wspace=0.3, hspace=0.6)
    ax1 = plt.subplot(grid[:10, 0])
    ax2 = plt.subplot(grid[:10, 1])
    ax3 = plt.subplot(grid[11, 0])
    ax4 = plt.subplot(grid[11, 1])

    sns.set(rc={'figure.figsize':(5,5)})
    sns.set_theme(style="white",palette=None)
    sns.set_style("white")

    sns.regplot(x=tmp_df['scores'].values.tolist(),
                y=tmp_df['exp'].values.tolist(),ax=ax1,color=color_mut)
    ax1.set_xlim(ax1.get_xlim()[0]-0.5,ax1.get_xlim()[1]+0.5)
    ax1.set_ylim(ax1.get_ylim()[0]-0.5,ax1.get_ylim()[1]+0.5)
    ax1.set_xlabel("PPI score")
    ax1.set_ylabel("mRNA expression")



    ax3.text(0.5 * (ax3.get_xlim()[0] + ax3.get_xlim()[1]),0,"R: "+str(np.round(r,4))+"   p-value: "+str(format(pval,'.1E')),fontsize=15,horizontalalignment='center',)

    scores_wt = all_wt_score_df[partner].values.tolist()
    gene_exp = df_wt_exp_samples.loc[gene].values.tolist()
    r = 0
    pval = 1
    try:
        if len(set(scores_wt))>3 and len(set(gene_exp))>3:
            outliers = get_otliers(scores_wt,gene_exp).index.values   
            tmp_df = pd.DataFrame()
            tmp_df['scores']=scores_wt
            tmp_df['exp']=gene_exp

            if len(outliers)>0:
                tmp_df.drop(tmp_df.index[outliers],axis=0,inplace=True)
            if len(set(tmp_df['scores'].values.tolist()))>2 and len(set(tmp_df['exp'].values.tolist()))>2:
                r,pval = ss.pearsonr(tmp_df['scores'].values.tolist(),tmp_df['exp'].values.tolist())
            N = len(tmp_df)
    except ValueError:
        r = 0
        pval = 1
        N = 0

    sns.regplot(x=tmp_df['scores'].values.tolist(),
                y=tmp_df['exp'].values.tolist(),ax=ax2,color=color_wt);
    ax2.set_xlim(ax2.get_xlim()[0]-0.5,ax2.get_xlim()[1]+0.5);
    ax2.set_ylim(ax2.get_ylim()[0]-0.5,ax2.get_ylim()[1]+0.5);
    ax2.set_xlabel("PPI score");
    ax2.set_ylabel("mRNA expression");
    ax4.text(0.5 * (ax4.get_xlim()[0] + ax4.get_xlim()[1]),0,"R: "+str(np.round(r,4))+"   p-value: "+str(format(pval,'.1E')),fontsize=15,horizontalalignment='center',)

    ax2.set_title(geneA+" WT/"+partner+ ": "+gene);
    ax1.set_title(geneA+" MUT/"+partner+ ": "+gene);

    for ax in [ax3,ax4]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])

    return(fig)




def pathway_enrichment(coding_genes,sig_genes,pathway_folder,files,savefile,partner):
    import statsmodels.sandbox.stats.multicomp as mltc
    import scipy.stats as stats
    import pandas as pd

    #all_genes = len(set(coding_genes['Approved symbol'].values.tolist()))
    all_genes = len(set(coding_genes))
    head = ["SET","fpval","num_sig_genes","num_sig_genes_included","num_genes_in_pathway",
                    "Expected, %","Actual, %","Enrichment","sig_genes_included"]
    enrichment_df = pd.DataFrame()
    bargraphs = []
    for geneset_file in files:
        pathway_f = open(pathway_folder+geneset_file, "r")
        #pathways = {}
        
        out = []
        for line in pathway_f:
            line = line.split("\t")
            
            #line = [x.strip().upper() for x in line[2:] if x.strip().upper() != '']
            #pathways[line[0]]=[x.strip().upper() for x in line[2:] if x.strip().upper() != '']
            name = line[0]
            pathway=[x.strip().upper() for x in line[2:] if x.strip().upper() != '']
            
    
            #for pathway in pathways.keys():
            #name = pathway
            #sig_genes_in_pathway = len(sig_genes.intersection(set(pathways[pathway])))
            
            sig_genes_in_pathway = len(sig_genes.intersection(set(pathway)))
            if sig_genes_in_pathway > 0:
                sig_genes_not_in_pathway = len(sig_genes)-sig_genes_in_pathway
        
                genes_in_pathways_not_sig_genes = len(pathway) - sig_genes_in_pathway
                genes_not_in_pathways_not_sig_genes = all_genes-sig_genes_in_pathway-sig_genes_not_in_pathway-genes_in_pathways_not_sig_genes
        
                exp = len(pathway)/all_genes*100
        
                if len(sig_genes)>0:
                    exp2 = sig_genes_in_pathway/len(sig_genes)*100
                else:
                    exp2 = 0
        
                ovr_status = ""
        
                m = [[sig_genes_in_pathway,genes_in_pathways_not_sig_genes],[sig_genes_not_in_pathway,genes_not_in_pathways_not_sig_genes]]

                oddratio, pval = stats.fisher_exact(m,alternative='greater')                              
        
                if exp == 0:
                    enrich = exp2
                elif exp2 == 0:
                    enrich = 0
                else:
                    enrich = exp2/exp
                nl = [name,pval,len(sig_genes),sig_genes_in_pathway,len(pathway),
                                exp,exp2,enrich,'; '.join(intersect(sig_genes,pathway))]
                out.append(nl)

        pathways_enrichment_df = pd.DataFrame(out,columns=head)    
        pvals = pathways_enrichment_df['fpval'].values
        if len(pvals)>0:
            qvals = mltc.multipletests(pvals,method='fdr_bh')[1]
        else:
            qvals = []
        pathways_enrichment_df['fqval']=qvals
        selected_df=pathways_enrichment_df.sort_values(by=['fqval'])
        bargraph = plt.figure(figsize=(3,3));
        plt.close(bargraph);
        if len(selected_df)>0:
            bargraph = enrichment_bargraph(selected_df, 10, geneset_file,savefile,partner);
            plt.close(bargraph)

        if len(enrichment_df)==0:
            enrichment_df = pathways_enrichment_df.copy()
        else:
            enrichment_df = pd.concat([enrichment_df,pathways_enrichment_df])
        bargraphs.append(bargraph)
        plt.close(bargraph)
    return enrichment_df,bargraphs


def enrichment_bargraph(selected_df, top, geneset_file,savefile,partner):
    import matplotlib.pyplot as plt
    import seaborn as sb

    #plt.close('all')
    plt.rcParams['pdf.fonttype'] = 42;
    plt.rcParams['ps.fonttype'] = 42;
    
    if top > len(selected_df):
        top = len(selected_df)
    selected_df.sort_values(by=['fqval'],inplace=True)
    fig, ax = plt.subplots(figsize=(3, (0.25*top)));
    
    
    grid = plt.GridSpec(12, 1, wspace=0.3, hspace=0.6);
    ax1 = plt.subplot(grid[0, 0]);
    ax2 = plt.subplot(grid[1:11, 0]);
    
    ax1.text(0,0,partner+" "+geneset_file+": Top-"+str(top)+" gene sets",fontsize=12,fontweight="bold");
    ax1.spines['top'].set_visible(False);
    ax1.spines['right'].set_visible(False);
    ax1.spines['left'].set_visible(False);
    ax1.spines['bottom'].set_visible(False);
    ax1.set_yticks([]);
    ax1.set_xticks([]);
    
    
    ax2.set_xlabel("-Log10(Q-Value)",fontsize=11);
    #ax2.set_title(geneset_file+": Top-"+str(top)+" gene sets",fontsize=15);
    sb.set_theme(style="white",palette=None);
    sb.set_style("whitegrid")
    plt.grid(False);
    mybar = sb.barplot(x=[-np.log10(x) for x in selected_df['fqval'].values[:top]],
               y = selected_df['SET'].values[:top],
               ax=ax2,dodge=False,palette='summer',orient='h');
    
    #mybar = sns.barplot(x=maths,y=physics,orient='h',ax=ax)
    mybar.set_yticklabels([]);
    ax2.set_yticks([]);
    lbl = selected_df['SET'].values[:top]
    mybar.bar_label(ax2.containers[0],label_type='edge',labels=lbl,fontsize=11,padding=5);
    ax2.spines['top'].set_visible(False);
    ax2.spines['right'].set_visible(False);
    
    #plt.vlines(x=-np.log10(0.05),
    #           ymin=0,
    #           ymax = ax2.get_ylim()[1],color='red')
    ax2.axvline(x=-np.log10(0.05),color='red',ls ='--',lw=0.5)
    #if savefile == 1:
    #    pdf_file = geneset_file+"_enrichment.png"
    #    fig.savefig(pdf_file);
            #, #dpi=600, #facecolor='w', edgecolor='w',
                                #orientation='portrait',
                                #format="pdf")#,
                                #bbox_inches="tight")#, 
                                #pad_inches=None,metadata=None)
    return fig;



def get_mut_data(mut_f,df_barcodes,cancer,geneA,mutA):
    mut_samples,wt_samples = get_mutation_data(mut_f, df_barcodes,cancer,geneA, mutA)
    mut_samples = list(set(mut_samples))
    wt_samples = list(set(wt_samples))
    #print("ALL_MUT_SAMPLES: ",len(mut_samples))
    #print("ALL_WT_SAMPLES: ",len(wt_samples))
    return(mut_samples,wt_samples)

def gen_filename(geneA,mutA,cancer,suffix):
    if len(mutA)==1:
        return(("_").join([geneA,mutA[0],cancer])+suffix)
    else: 
        return(("_").join([geneA,'MUT',cancer])+suffix)
    

def calculate_ppi_scores_not_scaled(out_folder,cancer,geneA,mutA,partners,df_mut_exp_samples,df_wt_exp_samples):
    from pathlib import Path
    import os
    #out_folder = out_folder+("_").join([cancer,geneA,("-").join(mutA),"batch"])
    #Path(out_folder).mkdir(parents=True, exist_ok=True)
    os.chdir(out_folder)


    all_mut_score_df = pd.DataFrame(index = df_mut_exp_samples.columns)
    all_wt_score_df = pd.DataFrame(index = df_wt_exp_samples.columns)
    for partner in partners:
        scores_mut,scores_mut_df = get_PPI_scores3_not_scaled(df_mut_exp_samples,geneA,partner)
        scores_wt,scores_wt_df = get_PPI_scores3_not_scaled(df_wt_exp_samples,geneA,partner)
        scores_mut_df.columns = [partner]
        scores_wt_df.columns = [partner]
        all_mut_score_df = pd.concat([all_mut_score_df,scores_mut_df],join='inner',axis=1)
        all_wt_score_df = pd.concat([all_wt_score_df,scores_wt_df],join='inner',axis=1)
    scores_mut_df.to_csv(out_folder+gen_filename(geneA,mutA,cancer,"_PPIscores.csv"),sep=",")
    scores_wt_df.to_csv(out_folder+geneA+"_WT_"+cancer+"_PPIscores.csv",sep=",")
    return(all_mut_score_df,all_wt_score_df,scores_mut_df,scores_wt_df)    
    

def get_coding_genes(hgnc_df):

    return(hgnc_df.loc[hgnc_df['Locus group']=='protein-coding gene'])
    
    

#Calculate PPI scores


def get_PPI_scores3_not_scaled(df1,gene1,gene2):
    exp1 = df1.loc[gene1]#.values
    exp2 = df1.loc[gene2]#.values
    exp1 = np.array(exp1.values.tolist())
    exp2 = np.array(exp2.values.tolist())
    
    
    
    scores = []
    for i,j in zip(exp1,exp2):
        scores.append(np.average([i,j]))
    samples = df1.columns.values.tolist()
    out_df = pd.DataFrame(index=samples)
    out_df['PPI_score']=scores
    return scores,out_df


  

#Functions for Enrichment analysis
#OVERREPRESENTATION oF SIGNATURE GENE SETS in HALLMARKS and PATHWAYS


def intersect(list1,list2):
    list1 = list(set(list1))
    list2 = list(set(list2))
    inter = [x for x in list1 if x in list2 and x.strip() != ""]
    return inter


def get_outliers(x,y):
    model = sm.OLS(x,y).fit()#formula=f, data=arr).fit()
    threshold = 4/len(x)
    influence = model.get_influence()
    #cooks = influence.cooks_distance
    (distance, p_value) = influence.cooks_distance
    influencial_data = distance[distance > threshold]
    return(influencial_data)


def get_otliers(x,y):
    tmp_df = pd.DataFrame()
    tmp_df['x']=x
    tmp_df['y']=y
    # fit the regression model using statsmodels library 
    f = 'x ~ y'
    model = ols(formula=f, data=tmp_df).fit()
    # calculate the cooks_distance - the OLSInfluence object contains multiple influence measurements
    cook_distance = OLSInfluence(model).cooks_distance
    (distance, p_value) = cook_distance
    threshold = 4/len(x)
    # the observations with Cook's distances higher than the threshold value are labeled in the plot
    influencial_data = distance[distance > threshold]
    return(influencial_data)

def compare_correlations(r1,r2,n1,n2):
    z1 = np.arctanh(r1)
    z2 = np.arctanh(r2)
    Z = (z1-z2) / (np.sqrt((1/(n1-3))+(1/(n2-3))))
    pval = ss.norm.sf(abs(Z))
    return pval
#get mutation data for the gene

def get_mutation_data(mut_f, df_barcodes,cancer,geneA, mutA):

    sample_IDs = df_barcodes.loc[:,cancer].values.tolist()
    sample_IDs = [x for x in sample_IDs if str(x)!='nan']
    mut_data = []
    head = ""
    N = 0
    with open(mut_f) as f:
        for line in f:
            line = line.split('\t')
            if N == 0:
                head = line
            if line[0].strip().upper() == geneA:
                mut_data.append(line)
            N +=1
    tmp = df_barcodes[cancer].dropna()
    print("Total samples:", len(tmp))
    df_gene_mut = pd.DataFrame(mut_data,columns = head)
    
    #samples with any mutation in geneA:
    samples_with_geneA_muts = df_gene_mut.loc[:,'Tumor_Sample_Barcode'].values.tolist()
    samples_with_geneA_muts = [x[:12] for x in samples_with_geneA_muts]
    
    if mutA != 'ALL':
        df_gene_mut = df_gene_mut.loc[df_gene_mut['HGVSp_Short'].isin(mutA)]
    mut_samples_all = df_gene_mut.loc[:,'Tumor_Sample_Barcode'].values.tolist() #all samples with the mutA mutation

    for i in range(len(mut_samples_all)):
        s = mut_samples_all[i].split("-")
        s = s[0]+"-"+s[1]+"-"+s[2]
        mut_samples_all[i] = s

    mut_samples = [] #mutated samples which also have expression data (included in sample_IDs list)

    for s in mut_samples_all:
        if s in sample_IDs:
            mut_samples.append(s)

    wt_samples = [] #wild type samples which also have expression data (included in sample_IDs list)
    for s in sample_IDs:
        #if s not in mut_samples_all and str(s) != 'nan':
        if s not in samples_with_geneA_muts and str(s) != 'nan':
            wt_samples.append(s)
    return(mut_samples,wt_samples)

def get_expression(exp_f,samples,coding_genes):
    df_exp = pd.read_csv(exp_f,sep='\t',index_col=0,header=0)
    new_columns = [i[:12] for i in df_exp.columns.values.tolist()]
    df_exp.columns = new_columns
    df_exp = df_exp.loc[:, ~df_exp.columns.duplicated()]
    #Keep samples for the cancer type only
    df_exp.drop([x for x in df_exp.columns if x not in samples],axis=1,inplace=True)
    #Keep protein coding genes only
    todrop = []
    df_exp.drop([x for x in df_exp.index if x not in coding_genes],inplace=True)
    
    #df_exp = keep_coding(df_exp,coding_genes)

    #Remove duplicated genes if any
    #df_exp.reset_index(inplace=True)
    todrop = []
    indices = []
    #for i in df_exp.index:
    #    if df_exp.iloc[i]['GENE_SYMBOL'] in indices:
    #        todrop.append(i)
    #    else:
    #        indices.append(df_exp.iloc[i]['GENE_SYMBOL'])
    for i in df_exp.index:
        if i in indices:
            todrop.append(i)
        else:
            indices.append(i)
    df_exp.drop(todrop,inplace=True,axis=0)
    #df_exp.set_index("GENE_SYMBOL",inplace=True)
    return df_exp

def drop_no_exp(df_exp):
    #print('Drop genes with zero expression in > 30% samples')
    todrop = []
    for c in df_exp.columns:
        zeros = len(df_exp.loc[df_exp[c]==0])
        if zeros*100/len(df_exp)>=30:
            todrop.append(c)
    df_exp.drop(todrop,axis=1,inplace=True)
    return df_exp

#Utilities
def prepare_cancer_barcodes(clinical_f):
    clinical_df = pd.read_csv(clinical_f,sep='\t')
    cancers = clinical_df['type'].unique().tolist()
    dict1 = {}
    for c in cancers:
        samples = clinical_df.loc[clinical_df['type']==c]['bcr_patient_barcode'].values.tolist()
        dict1[c]=samples
    barcode_df = pd.DataFrame.from_dict(dict1,orient='index').transpose()
    return barcode_df


def get_wt_mut_expression(mut_f,exp_f,coding_genes,barcode_df,cancer,geneA,mutA,out_folder):
    mut_samples,wt_samples = get_mut_data(mut_f,barcode_df,cancer,geneA,mutA)
    df_exp = get_expression(exp_f,mut_samples+wt_samples, coding_genes)

    #Drop genes with zero expression in > 30% samples in up_df
    df_exp = drop_no_exp(df_exp.copy())#av

    #Log(x+1) taransform
    df_exp = df_exp.apply(lambda x: np.log2(x+1),)
    
    df_wt_exp_samples = df_exp.drop([x for x in mut_samples if x in df_exp.columns],axis=1)
    df_mut_exp_samples = df_exp.drop([x for x in wt_samples if x in df_exp.columns],axis=1)

    df_wt_exp_samples.to_csv(out_folder+geneA+"_WT_"+cancer+"_exp.csv",sep=",")
    df_mut_exp_samples.to_csv(out_folder+gen_filename(geneA,mutA,cancer,"_exp.csv"))
    return(df_mut_exp_samples,df_wt_exp_samples)


def get_mRNA_expession(exp_f):
    
    genes_with_expression = []
    n = 0
    with open(exp_f) as f:
        for line in f:
            line = line.strip()
            if line != '':
                line = line.replace('"','').split('\t')
            if n > 0 and line != '':
                symbol = line[0]
                genes_with_expression.append(symbol)
            n = 1
    genes_with_expression = list(set(genes_with_expression))
    return(genes_with_expression)


def get_cancer_mutation_freq(barcode_df,mut_f,geneA,mutA,cancers):
#Mutation frequency
    n = 0
    cancers_df = pd.DataFrame(index = cancers,columns = ['Wild type','Other mutants','Target mutant'])
    cancers_df['Target mutant']=0
    cancers_df['Other mutants']=0
    for c in cancers_df.index:
        cancers_df['Wild type'].loc[c] = len(barcode_df[c].dropna())
    with open(mut_f) as f:
        used_codes = []
        for line in f:
            line = line.split("\t")
            if n == 0:
                ind = line.index('Tumor_Sample_Barcode') #the column number for the Sample ID
                n=1
            if geneA in line:
                flag2 = 0 #indicate the presence of the target mutatioin(s) 
                for m in mutA:
                    if m in line:
                        flag2 = 1
                        break
                barcode = line[ind][:12]
                if barcode not in used_codes:
                    used_codes.append(barcode)
                    for c in barcode_df.columns:
                        if barcode in barcode_df[c].values:
                            if flag2 == 1:
                                cancers_df['Target mutant'].loc[c] = cancers_df['Target mutant'].loc[c]+1
                            else:
                                cancers_df['Other mutants'].loc[c] = cancers_df['Other mutants'].loc[c]+1
                            cancers_df['Wild type'].loc[c] = cancers_df['Wild type'].loc[c]-1
        
    cancers_df.sort_values(by=["Target mutant"],ascending=False,inplace=True)
    fig, ax = plt.subplots(figsize=(10, 4))
    cancers_df.plot(kind='bar', stacked=True, color=['lightgrey', 'skyblue', 'green'],ax=ax,
                   width=0.8)
    ax.set_ylabel("Number of samples")
    ax.set_xlabel("Cancer type")
    # Add labels to each bar.
    y_offset = 20
    for i, target_mut in enumerate(cancers_df['Target mutant'].values.tolist()):
        total = sum(cancers_df.iloc[i])
        cancer = cancers_df.index[i]
        if target_mut > 0:
            ax.text(i, total + y_offset, round(target_mut), ha='center')
    return (fig, cancers_df)


def get_binding_partners(ppi_f):
    partners = []
    with open(ppi_f) as f:
        for partner in f:
            partner = partner.strip()
            partners.append(partner)
    return(list(set(partners)))

def compare_ppi_scores(ppi_score_dict,partners):
    out = []
    header = []
    for i in range(0,len(list(ppi_score_dict.keys()))):
        c1 = list(ppi_score_dict.keys())[i]
        header.append(c1+"_AVR")
        header.append(c1+"_SD")
        for j in range(i+1,len(list(ppi_score_dict.keys()))):
            c2 = list(ppi_score_dict.keys())[j]
            cc = "PVAL_"+c1+"-"+c2
            cc_fc = "FC_"+c1+"-"+c2
            header.append(cc)
            header.append(cc_fc)
    header.append("KW_H_Test_PVAL")
    out = []
    for gene in partners:
        ppi_scores = []
        for cancer in ppi_score_dict.keys():
            ppi_scores.append(ppi_score_dict[cancer][0][gene].values)
    
    #Kruskal-Wallis H-test
        kstat, pval = ss.kruskal(*ppi_scores)
        p_values= sp.posthoc_dunn(ppi_scores, p_adjust = 'holm')
        p_values.columns = ppi_score_dict.keys()
        p_values['cancer'] = list(ppi_score_dict.keys())
        p_values.set_index('cancer',inplace=True)
        nl = []
        for i in range(0,len(list(ppi_score_dict.keys()))):
            c1 = list(ppi_score_dict.keys())[i]
            avr1 = np.mean(ppi_score_dict[c1][0][gene])
            sd1 = np.std(ppi_score_dict[c1][0][gene])
            nl.append(avr1)
            nl.append(sd1)
            for j in range(i+1,len(list(ppi_score_dict.keys()))):
                c2 = list(ppi_score_dict.keys())[j]
                cc = "PVAL_"+c1+"-"+c2
                cc_value = p_values[c1].loc[c2]
                avr2 = np.mean(ppi_score_dict[c2][0][gene])
                sd2 = np.std(ppi_score_dict[c2][0][gene])
                cc_fc = "FC_"+c1+"-"+c2
                fc_value = avr1/avr2
                nl.append(cc_value)
                nl.append(fc_value)
        nl.append(pval)
        out.append(nl)
    out_df = pd.DataFrame(out,columns=header)       
    out_df['neoPPI'] = partners
    out_df.set_index('neoPPI',inplace=True)
    #Add FDR
    out_df1 = pd.DataFrame(columns = out_df.columns) #KW_H_Test_PVAL != NaN
    out_df2 = pd.DataFrame(columns = out_df.columns) #KW_H_Test_PVAL == NaN
    for i in out_df.index:
        if str(out_df.loc[i]['KW_H_Test_PVAL']) =='nan':
            out_df2.loc[i] = out_df.loc[i]
        else:
            out_df1.loc[i] = out_df.loc[i]
    out_df1['KW_H_Test_QVAL'] = mltc.multipletests(out_df1['KW_H_Test_PVAL'].values.tolist(),method='fdr_bh')[1]
    out_df2['KW_H_Test_QVAL'] = np.nan
    
    out_df = pd.concat([out_df1,out_df2],join='inner')
    return(out_df)

def get_neo_ppi_score_heatmap(out_df):
    map_df = out_df.copy()
    map_df.drop([x for x in map_df.columns if "AVR" not in x],axis=1,inplace=True)
    map_df.columns = [x.split("_")[0] for x in map_df.columns]
    map_df = map_df.transpose()
    g = sns.clustermap(data = map_df,figsize=(20,4),dendrogram_ratio=(.025, 0.2),xticklabels=1,cbar_pos=(-0.02, 0.75, .01, .2),
                      cmap='viridis',method='ward')
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 8);
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 10,va="center");
    return(g)

def build_heatmap(cancer,ppi_score_dict,clinical_f,geneA,mutA,fig_folder,sortby):
    
    df = ppi_score_dict[cancer][0].copy()
    df.fillna(0,inplace=True)
    gene_exp = ppi_score_dict[cancer][1].loc[geneA].values
    gene_tmp1 = gene_exp.copy()
    gene_tmp1.sort()
    colors = plt.get_cmap('viridis')(np.linspace(0.1, 0.9, len(gene_tmp1)))
    lut = dict()
    for i in range(0,len(gene_tmp1)):
        lut[gene_tmp1[i]]=colors[i]
    
    gene_df = ppi_score_dict[cancer][1].transpose().pop(geneA)
    race_df = pd.read_csv(clinical_f,sep="\t",index_col=0).loc[df.index].pop("race")
    gender_df = pd.read_csv(clinical_f,sep="\t",index_col=0).loc[df.index].pop("gender")
    age_df = pd.read_csv(clinical_f,sep="\t",index_col=0).loc[df.index].pop("age_at_initial_pathologic_diagnosis")
    age_df.fillna(0,inplace=True)

    race_lut = dict()
    for i in race_df.values:
        if i =="WHITE":
            race_lut[i]="whitesmoke"
        elif i == "BLACK":
            race_lut[i]="BLACK"
        elif i == "ASIAN":
            race_lut[i] = "Orange"
        else:
            race_lut[i] = "cyan" 
        
    gender_lut = dict(zip(gender_df.unique(), ["deepskyblue","lightpink"]))

    age = age_df.values.tolist()
    #age_lut = dict()
    #for i in range(0,len(age)):
    #    if age[i] < 40:
    #        c = 'red'
    #    elif age[i] >= 40 and age[i] < 50:
    #        c = 'orange'
    #    else:
    #        c = 'green'
    age_pal = sns.light_palette('skyblue', len(set(age)))
    age_lut = dict(zip(set(age), age_pal))

    gene_mut_name=geneA+"_"+mutA[0]
    row_colors = pd.DataFrame({#gene_mut_name:gene_df.map(lut),
                           'Race':race_df.map(race_lut),
                           'Gender':gender_df.map(gender_lut),
                          'Age':age_df.map(age_lut)})

    tmp1 = pd.DataFrame({#gene_mut_name:gene_df,
                           'Race':race_df,
                           'Gender':gender_df,
                          'Age':age_df})
    
    df = pd.concat([row_colors,df],axis=1,join='inner')
    df.sort_values(by=sortby,inplace=True,ascending=False)
    df.drop(['Age','Race','Gender'],inplace=True,axis=1)
    
    df1 = pd.concat([tmp1,df],axis=1,join='inner')
    df1.sort_values(by=sortby,inplace=True,ascending=False)
    
    max_val = max(df.max().values)
    min_val = min(df.min().values)
    
    
    if len(df.columns)>1:
        g = sns.clustermap(data=df,cmap='seismic',method='complete',row_colors=row_colors,xticklabels=1,yticklabels=0,
                  dendrogram_ratio=(.1, 0.1),row_cluster=False,figsize=(6,6),
                       cbar_kws=dict(ticks=[min_val, max_val], orientation='horizontal'))
        g.ax_col_dendrogram.set_title(cancer)
    else:
        g = sns.clustermap(data=df,cmap='seismic',method='complete',row_colors=row_colors,xticklabels=1,yticklabels=0,
                  dendrogram_ratio=(.1, 0.1),row_cluster=False,col_cluster=False,figsize=(6,6),
                       cbar_kws=dict(ticks=[min_val, max_val], orientation='horizontal'))
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 6);
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 3);
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    #g.fig.suptitle(cancer)
    
    x0, _y0, _w, _h = g.cbar_pos
    g.ax_cbar.set_position([x0, 0.95, g.ax_row_dendrogram.get_position().width, 0.02])
    fig_file = fig_folder+gen_filename(geneA,mutA,cancer,"heatmap.pdf")
    g.savefig(fig_file, dpi=600,format='pdf')
    return(g,df1)

def get_best_grid(total,columns):
    rem0 = total
    rows_to_use = 1
    cols_to_use = 1
    while columns > 1:
        rem = total % columns
        rows = int(total/columns)
        
        if rem == 0:
            rows_to_use = rows
            cols_to_use = columns
            break
        else:
            if rem < rem0:
                rem0 = rem
                rows_to_use = rows
                cols_to_use = columns
        columns-=1
    return rows_to_use, cols_to_use


def boxplot_neoppi_score_distribution(genes,cancers,ppi_score_dict,palette,color):
    sns.set_style("whitegrid", {'axes.grid' : False})
    rows,cols = get_best_grid(len(genes),5)
    plot_width = 3.4
    plot_height = 3.2
    fig_size = (cols*plot_width, rows*plot_height)
    fig,axs = plt.subplots(rows,cols,figsize=fig_size)

    fig.tight_layout(h_pad=2)
    i = 0
    j = 0
    for gene in genes:
        gene_df = pd.DataFrame(columns = ['PPIscore','cancer'])
        for cancer in cancers:
            cdf = pd.DataFrame(ppi_score_dict[cancer][0][gene].values.tolist(),columns=['PPIscore'])
            cdf['cancer']=cancer
            gene_df = pd.concat([gene_df,cdf],axis=0)
        gene_df.reset_index(inplace=True)  
    
        if rows > 1 and cols > 1:
            ax = axs[i,j]
        elif rows == 1 and cols > 1:
            ax = axs[j]
        elif rows > 1 and cols == 1:
            ax = axs[i]
        else:
            ax = axs
        ax.set_title(gene)
        
        #sns.swarmplot(data = gene_df,ax=ax,size=2)
        #sns.boxplot(data = gene_df,ax=ax,palette='Set3',fliersize=0,linewidth=0.5)
        sns.swarmplot(data = gene_df,x = 'cancer', y = 'PPIscore',size=2,color=color,ax=ax)
        sns.boxplot(data = gene_df,x = 'cancer', y = 'PPIscore',palette=palette,fliersize=0,linewidth=0.5,ax=ax)
        ax.set_ylabel("PPI score")
        ax.set(xlabel=None)
        ax.set_yticks(np.arange(int(ax.get_ylim()[0]), int(ax.get_ylim()[1])+1, 2))
        
        j+=1
        if j == cols:
            j=0
            i+=1
    return(fig)
        

def survival_analysis_with_progress_bar(mut_samples,clinical_f,all_mut_score_nn_df,show_plot,pval_treshold):
    '''
    Determine the difference in survival times in mutant samples with high and low neoPPI scores
    mut_samples: list 
                 A list of patient TCGA IDs for samples with the mutation
    clinical_f: str
                path to the file with clinical data
    all_mut_score_nn_df: 
                DataFrame
                A DataFrame with neoPPI scores
    show_plot: str
               Indicates if the survival plots should be shown for 'all', 'none', or 'significant' correlations
    pval_treshold: float
                The p-value treshold fo statistical significance
    
    '''
    if show_plot not in ['all','significant','none']:
        raise Exception("show_plot parameter should be either 'all', 'significant', or 'none'")
        exit()

    mut_survive_df = get_mut_survival_data(clinical_f,mut_samples)[['bcr_patient_barcode','OS.time','vital_status']]
    mut_survive_df['vital_status'].replace(['Alive', 'Dead'],[0, 1], inplace=True)
    mut_survive_df.set_index("bcr_patient_barcode",inplace=True)
    mut_survive_df = pd.concat([mut_survive_df,all_mut_score_nn_df],join='inner',axis=1)
    mut_survive_df['OS.time'] = [(float(x)*0.0328767)/12 for x in mut_survive_df['OS.time'].values.tolist()]

    out = []
    fig_dict = {}
    max_count = len(mut_survive_df.columns[2:])
    progressbar = IntProgress(min=0, max=max_count)
    display(progressbar)
    
    for partner in mut_survive_df.columns[2:]:
        partner_df = mut_survive_df[['OS.time','vital_status',partner]].copy()
        partner_df.dropna(inplace=True)
        T = partner_df["OS.time"]
        E = partner_df["vital_status"]
        kmf1 = KaplanMeierFitter()
        kmf1.fit(durations = T, event_observed = E)
        kmf2 = KaplanMeierFitter()
        kmf2.fit(durations = T, event_observed = E)
    
        m1 = (partner_df[partner] >= np.median(partner_df[partner]))
        m2 = (partner_df[partner] < np.median(partner_df[partner]))
        pval = 1
        try:
            results = logrank_test(T[m1], T[m2], E[m1], E[m2], alpha=.99)
            pval = results.p_value
            sns.set_theme(style="white",palette=None)
            sns.set_style("white")
        except:
            pass
        kmf1.fit(durations = T[m1], event_observed = E[m1], label = "High")
        mtime1 = median_survival_times(kmf1.survival_function_)
        kmf2.fit(durations = T[m2], event_observed = E[m2], label = "Low")
        mtime2 = median_survival_times(kmf2.survival_function_)
        out.append([partner,mtime1,mtime2,pval])
        flag = 0
        if show_plot == 'significant':
            if (mtime1 < mtime2) and (pval < pval_treshold):
                flag = 1
            else:
                flag = 0
        elif show_plot == 'none':
            flag = 0
        else:
            flag = 1
        if flag == 1:
            fig,ax = plt.subplots(figsize=(4,4))
            kmf1.survival_function_.plot(ax=ax,color='red')
            kmf2.survival_function_.plot(ax=ax)
            ax.set_title(partner)
            ax.set_xlabel("Years")
            ax.set_ylabel("Survival probability")
            ax.text(0,0.1,"P-value = "+str(np.round(pval,3)))
            fig_dict[partner] = [kmf1,kmf2,pval]
            plt.close(fig)
        
        progressbar.value += 1
    surv_sum_df = pd.DataFrame(out,columns = ['PARTNER','MEDIAN_TIME_HIGH','MEDIAN_TIME_LOW','PVALUE'])
    
    plots = len(fig_dict.keys())
    if plots > 0:
        rows,cols = get_best_grid(plots,5)
        if plots/cols > rows:
            rows = rows + 1
        fig, axs = plt.subplots(rows,cols,figsize=(cols*4,rows*3.5))
        fig.tight_layout(pad=4)
        i = 0
        j = 0
        for partner in fig_dict.keys():
            if rows > 1 and cols > 1:
                ax = axs[i,j]
            elif rows == 1 and cols > 1:
                ax = axs[j]
            elif rows > 1 and cols == 1:
                ax = axs[i]
            else:
                ax = axs
            #fig_dict[partner][0].survival_function_.plot(ax=axs[i,j],color='red')
            #fig_dict[partner][1].survival_function_.plot(ax=axs[i,j])
            fig_dict[partner][0].survival_function_.plot(ax=ax,color='red')
            fig_dict[partner][1].survival_function_.plot(ax=ax)
        
            ax.set_title(partner)
            ax.set_xlabel("Years")
            ax.set_ylabel("Survival probability")
            ax.set_ylim(0,1)
            ax.text(0,0.05,"P-value = "+str(np.round(fig_dict[partner][2],3)),fontsize=8)
            j+=1
            if j == cols:
                j = 0
                i = i + 1
    else:
        fig = plt.figure();
        plt.close(fig)
    return(surv_sum_df,fig)

def calculate_correlations(df_mut_exp_samples,df_wt_exp_samples,partners,all_mut_score_nn_df,all_wt_score_nn_df):
    exp_arr_mut = df_mut_exp_samples.values
    exp_genes_mut = df_mut_exp_samples.index

    exp_arr_wt = df_wt_exp_samples.values
    exp_genes_wt = df_wt_exp_samples.index

    corr_dict = {}
    
    max_count = len(partners)
    progressbar = IntProgress(min=0, max=max_count)
    display(progressbar)
    
    for partner in partners:
        scores_mut = all_mut_score_nn_df[partner].values

        out = get_corr_genes(scores_mut,exp_arr_mut,exp_genes_mut)
    
        out_mut_df = pd.DataFrame(out,columns = ['GENE','CORR_BP_MUT','PVAL_BP_MUT','N_MUT'])
        out_mut_df.set_index('GENE',inplace=True)
        #For the WT samples:
        scores_wt = all_wt_score_nn_df[partner].values
        out = get_corr_genes(scores_wt,exp_arr_wt,exp_genes_wt)
        out_wt_df = pd.DataFrame(out,columns = ['GENE','CORR_BP_WT','PVAL_BP_WT','N_WT'])
        out_wt_df.set_index('GENE',inplace=True)
        #Combine:
        out_df = pd.concat([out_mut_df,out_wt_df],axis=1)
        #Append to the dictionary
        corr_dict[partner] = out_df
        progressbar.value += 1
    #Get FDR
    for partner in corr_dict.keys():
        pval_mut_df = pd.DataFrame(corr_dict[partner]['PVAL_BP_MUT'])
        pval_mut_df.dropna(inplace=True)
        pval_mut_df['QVAL_BP_MUT'] = mltc.multipletests(pval_mut_df['PVAL_BP_MUT'].values.tolist(),method='fdr_bh')[1]
    
        pval_wt_df = pd.DataFrame(corr_dict[partner]['PVAL_BP_WT'])
        pval_wt_df.dropna(inplace=True)
        pval_wt_df['QVAL_BP_WT'] = mltc.multipletests(pval_wt_df['PVAL_BP_WT'].values.tolist(),method='fdr_bh')[1]
    
        corr_dict[partner].drop(['PVAL_BP_MUT','PVAL_BP_WT'],inplace=True,axis=1)
    
        corr_dict[partner] = pd.concat([corr_dict[partner],pval_mut_df,pval_wt_df],axis=1)
    
        #Compare the MUT and WT correlations
        corr_dict[partner]['PVAL'] = compare_correlations(corr_dict[partner]['CORR_BP_MUT'],
                                                         corr_dict[partner]['CORR_BP_WT'],
                                                         corr_dict[partner]['N_MUT'],
                                                        corr_dict[partner]['N_WT'])
        pval_df = pd.DataFrame(corr_dict[partner]['PVAL'])
        pval_df.dropna(inplace=True)
        pval_df['QVAL'] = mltc.multipletests(pval_df['PVAL'].values.tolist(),method='fdr_bh')[1]
    
        corr_dict[partner].drop(['PVAL'],inplace=True,axis=1)
        corr_dict[partner] = pd.concat([corr_dict[partner],pval_df],axis=1)
    return(corr_dict)
    
sig_genes = []

def on_value_change(change):
    global sig_genes
    sig_genes, tbl = get_signature_genes2(corr_dict,partner,geneA,box.children[0].children[0].value,
                    box.children[0].children[1].value,box.children[0].children[2].value,box.children[0].children[3].value)
    return(box,sig_genes)

def get_signature_genes_for_multiple_binding_partners(CORR_BP_MUT,PVAL_BP_MUT,PVAL,QVAL,corr_dict,cancer,geneA,mutA):
    sign_gene_dict = {}
    for p in corr_dict.keys():
        sign_gene_dict[p] = corr_dict[p].loc[(corr_dict[p]['CORR_BP_MUT']>=CORR_BP_MUT) &
                                            (corr_dict[p]['PVAL_BP_MUT']<=PVAL_BP_MUT) &
                                            (corr_dict[p]['CORR_BP_MUT']>corr_dict[p]['CORR_BP_WT']) &
                                            (corr_dict[p]['PVAL']<=PVAL) &
                                            (corr_dict[p]['QVAL']<=QVAL)]
        sign_gene_dict[p]['CANCER']=cancer
        sign_gene_dict[p]['DRIVER']=geneA+"_"+("_").join(mutA)
        sign_gene_dict[p]['PARTNER']=p
    sign_gene_df = pd.DataFrame(columns = sign_gene_dict[list(sign_gene_dict.keys())[0]].columns)
    for p in sign_gene_dict.keys():
        sign_gene_df = pd.concat([sign_gene_df,sign_gene_dict[p]],axis=0)
    return(sign_gene_dict,sign_gene_df)

def sig_genes_size(sign_gene_dict):
    out = []
    for p in sign_gene_dict.keys():
        out.append([p,len(sign_gene_dict[p].index)])
    df_stat= pd.DataFrame(out,columns=["neoPPI","SIZE"])
    df_stat.set_index("neoPPI",inplace=True)
    df_stat.sort_index(inplace=True)
    return(df_stat)
    
def jaccard(sign_gene_dict,fig_folder):
    root = tkinter.Tk();
    genesets = []
    genesets_data = []
    all_J_scores = []
    jout = []
    Jarr = []

    sign_gene_dict1 = {}
    for i in sign_gene_dict.keys():
        if len(sign_gene_dict[i].index)>0:
            sign_gene_dict1[i] = sign_gene_dict[i]

    for i in sign_gene_dict1.keys():#range(len(partners)):
        partner1 = i
        genes1 = sign_gene_dict1[i].index.tolist()
        if len(genes1)>0:
            s = len(genes1)
            kk = [partner1,s]
            genesets.append(kk)
            genesets_data.append(s)
            row_J_scores = []
            for j in sign_gene_dict1.keys():
                partner2 = j
                genes2 = sign_gene_dict1[j].index.tolist()
                T = len(genes1)
                C = 0
                for g2 in genes2:
                    if g2 not in genes1:
                        T = T + 1
                    else:
                        C = C +1
                J = C/T*100
                if i != j:
                    Jarr.append(J)
                row_J_scores.append(J)
            all_J_scores.append(row_J_scores)
    jdf = pd.DataFrame(all_J_scores,index = sign_gene_dict1.keys(), columns = sign_gene_dict1.keys())
    
    if len(jdf.columns) < 2:
        print("Cannot build a cluster map for one partner")
    else:
    #adjust heatmap size:
        coeff = 1/root.winfo_fpixels('1i')
        width = len(jdf.columns.tolist())*coeff*100
        height = len(jdf.index.tolist())*coeff*100
        if width > 2^16 or height > 2^16:
            width = 2^16
            height = 2^16
    g = sns.clustermap(jdf,figsize=(width,height),
                    xticklabels=1,yticklabels=1,
                    cbar_pos=(-0.04, 0.75, .01, .2),
                    dendrogram_ratio=(.025, 0.025),
                    cmap='Blues')
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 4);
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 4);
    print("Average Jaccard Index: ", np.round(np.mean(Jarr),1), "%")
    #g.savefig(fig_folder+"jackard.png",dpi=600)
    return(g)


def neoppi_sig_gene_heatmap(sign_gene_df):
    root=tkinter.Tk()
    sub_sign_gene_df = sign_gene_df[['CORR_BP_MUT','PARTNER']].copy()
    sub_sign_gene_df.reset_index(inplace=True)
    sub_sign_gene_df.columns=['SIG_GENE','CORR_BP_MUT','PARTNER']
    sig_gene = pd.DataFrame(index = list(set(sub_sign_gene_df['PARTNER'])),columns = list(set(sub_sign_gene_df['SIG_GENE'])))
    for i in sub_sign_gene_df.index:
        g = sub_sign_gene_df.loc[i]['SIG_GENE']
        p = sub_sign_gene_df.loc[i]['PARTNER']
        sig_gene[g].loc[p]=sub_sign_gene_df.loc[i]['CORR_BP_MUT']
    sig_gene.fillna(0,inplace=True)
    
    if len(sig_gene.columns) < 2:
        print("Cannot build a cluster map for one partner")
    else:
    #adjust heatmap size:
        coeff = 1/root.winfo_fpixels('1i')
        width = len(sig_gene.columns.tolist())*coeff*20
        height = len(sig_gene.index.tolist())*coeff*500
        if width > 2^16 or height > 2^16:
            width = 2^16
            height = 2^16
    
    g = sns.clustermap(sig_gene,figsize=(width,height),
                    xticklabels=1,yticklabels=1,
                    cbar_pos=(-0.04, 0.75, .01, .2),
                    dendrogram_ratio=(.025, 0.025),
                    cmap='Blues',method='ward')
    if len(sig_gene.columns.tolist())>50:
        fontsize = 5
    elif len(sig_gene.columns.tolist())>20:
        fontsize = 8
    else:
        fontsize = 10
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = fontsize);
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize =fontsize);
    g.ax_heatmap.set_xlabel("Signature gene");
    g.ax_heatmap.set_ylabel("neoPPI");
    
    return(g)

def pathway(sign_gene_dict,files,partner,coding_genes,pathway_folder):
    partner_genes = []
    bars_dict = {}
    enrichment_dict = {}
    for partner in list(sign_gene_dict.keys()):
        sig_genes = set(sign_gene_dict[partner].index.tolist())
        enrichment_df,bargraph = pathway_enrichment(coding_genes,sig_genes,pathway_folder,files,1,partner);
        enrichment_df.reset_index(inplace=True,drop=True)
        bars_dict[partner]=dict();#bargraph;
        for i in range(0,len(files)):
            bars_dict[partner][files[i]]=bargraph[i]
        enrichment_dict[partner] = enrichment_df
        enriched_pathways = enrichment_df.loc[(enrichment_df['fqval']<0.05)]
        if len(enriched_pathways) > 0:
            partner_genes.append(partner)
        for i in enriched_pathways.index:
            pset = enriched_pathways.loc[i]['SET']
        partner_genes = list(set(partner_genes))
    return(bars_dict,enrichment_dict)

def connect_mutant_driver_partner_pathway(enrichment_dict,geneA):
    net_arr = []
    for p in enrichment_dict.keys():
        tmp_df = enrichment_dict[p].copy()
        tmp_df = tmp_df.loc[tmp_df['fqval']<0.05]
        net_arr.append([p,geneA+"_MUT"])
        for i in tmp_df.index:
            net_arr.append([p,tmp_df.loc[i]['SET']])
    all_enrichment_df = pd.DataFrame(net_arr,columns=["Partner","Pathway"])
    net_types_arr = []
    partner_genes = all_enrichment_df['Partner'].unique().tolist()
    pathways = [x for x in all_enrichment_df['Pathway'].unique().tolist() if x != geneA+"_MUT"]

    for i in partner_genes:
        net_types_arr.append([i,"Partner"])
    for i in pathways:
        net_types_arr.append([i,"Pathway"])
    net_types_arr.append([geneA+"_MUT","Driver"])
    
    all_enrichment_types_df = pd.DataFrame(net_types_arr,columns=['node','type'])
    return(all_enrichment_df,all_enrichment_types_df)

def enrichment_heatmap(subset,all_enrichment_df,geneA,fig_folder):
    root = tkinter.Tk();
    todrop = [x for x in all_enrichment_df.index if subset not in all_enrichment_df.loc[x]['Pathway'] and
            geneA+"_MUT" != all_enrichment_df.loc[x]['Partner']]
    subset_enrichment_df = all_enrichment_df.drop(todrop)

    enriched_partners = subset_enrichment_df['Partner'].unique().tolist()
    enriched_pathways = subset_enrichment_df['Pathway'].unique().tolist()
    enriched_map = pd.DataFrame(index = enriched_pathways,columns=enriched_partners)
    for c in enriched_map.columns:
        enriched_map[c]=0
    for p in enriched_map.index:
        pdf = subset_enrichment_df.loc[subset_enrichment_df['Pathway']==p]
        pgenes = pdf['Partner'].tolist()
        for g in pgenes:
            enriched_map[g].loc[p]=1
    if len(enriched_map.columns) < 2:
        print("Cannot build a cluster map for one partner")
    else:
    #adjust heatmap size:
        coeff = 1/root.winfo_fpixels('1i')
        width = len(enriched_map.columns.tolist())*coeff*20
        height = len(enriched_map.index.tolist())*coeff*12
    
    
        g = sns.clustermap(data=enriched_map,cmap=['white','green'],
                    xticklabels=1,yticklabels=1,method='ward',
                    dendrogram_ratio=(.05, 0.05),cbar_pos=None,
                    linewidths=0.5, linecolor='lightgrey',figsize=(width,height))
    
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 8);
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 8);
        #border around the map:
        for _, spine in g.ax_heatmap.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(1)
        
    return(g)
def single_gene_survival_plot(df_mut_exp_samples,gene,clinical_f):
    mut_samples = df_mut_exp_samples.columns
    out=[]
    if len(mut_samples) > 0:
        #sig_genes = set(sign_gene_dict[partner].index.tolist())
        mut_survive_df = get_mut_survival_data(clinical_f,mut_samples)
        gene = gene.upper()
        exp_cutoff=33
        pval=0.05
        savefile=0
        df_mut_high,df_mut_low,target_gene = get_high_low_samples(df_mut_exp_samples,gene,exp_cutoff)
        mut_high_samples = [x[0:12] for x in df_mut_high.columns.values.tolist()]
        mut_low_samples = [x[0:12] for x in df_mut_low.columns.values.tolist()]
        mut_high_surv_df = mut_survive_df.loc[mut_survive_df['bcr_patient_barcode'].isin(mut_high_samples)]
        mut_low_surv_df = mut_survive_df.loc[mut_survive_df['bcr_patient_barcode'].isin(mut_low_samples)]
        
        p_value,med_time_high,med_time_low,fig = survival_times("",gene,mut_high_surv_df,mut_low_surv_df,pval,savefile)
        out = [fig,p_value,med_time_high,med_time_low]
    return(out)

def get_sig_gene_survival_plot(partner,df_mut_exp_samples,sign_gene_dict,clinical_f):
    mut_samples = df_mut_exp_samples.columns
    figs = {}
    if len(mut_samples) > 0:
        sig_genes = set(sign_gene_dict[partner].index.tolist())
        mut_survive_df = get_mut_survival_data(clinical_f,mut_samples)
        surv_df,figs = conduct_survival_analysis(partner,sig_genes,df_mut_exp_samples,mut_survive_df,
                                                 exp_cutoff=33,pval=0.05,savefile=1)
        surv_df.set_index('GENE',inplace=True)
        cols = surv_df.columns.tolist()
        cols[-1] = 'CLIN_PVAL'
        surv_df.columns = cols
        #Add survival data to signature gene annotation
        sign_gene_dict[partner] = pd.concat([sign_gene_dict[partner],surv_df],axis=1)
        #print(len(surv_df[surv_df['CLIN_PVAL']<0.05]))
    return(sign_gene_dict,figs)

def conduct_survival_analysis(partner, sig_genes,df_mut_exp_samples,mut_survive_df,exp_cutoff,pval,savefile):
    out = []
    figs = {}
    for sig_gene in sig_genes:
        sig_gene = sig_gene.upper()
        df_mut_high,df_mut_low,target_gene = get_high_low_samples(df_mut_exp_samples,sig_gene,exp_cutoff)
        mut_high_samples = [x[0:12] for x in df_mut_high.columns.values.tolist()]
        mut_low_samples = [x[0:12] for x in df_mut_low.columns.values.tolist()]
        mut_high_surv_df = mut_survive_df.loc[mut_survive_df['bcr_patient_barcode'].isin(mut_high_samples)]
        mut_low_surv_df = mut_survive_df.loc[mut_survive_df['bcr_patient_barcode'].isin(mut_low_samples)]
        
        p_value,med_time_high,med_time_low,fig = survival_times(partner,sig_gene,mut_high_surv_df,mut_low_surv_df,pval,savefile)
    
        figs[sig_gene] = [fig,p_value,med_time_high,med_time_low]
        out.append([sig_gene,med_time_high,med_time_low,p_value])
    surv_df = pd.DataFrame(out,columns=['GENE','med_time_high','med_time_low','PVAL'])
    return(surv_df,figs)

def survival_times(partner,sig_gene,mut_survive_df,wt_survive_df,pval,savefile):
    #former survival
    gr1_status = mut_survive_df['vital_status'].values.tolist()
    gr1_status = [1 if x == 'Dead' else 0 for x in gr1_status]
    gr1_times = mut_survive_df['OS.time'].values.tolist()

    gr2_status = wt_survive_df['vital_status'].values.tolist()
    gr2_status = [1 if x == 'Dead' else 0 for x in gr2_status]
    gr2_times = wt_survive_df['OS.time'].values.tolist()

    #convert days to months
    gr1_times = [(float(x)*0.03287671)/12 for x in gr1_times]
    gr2_times = [(float(x)*0.03287671)/12 for x in gr2_times]
    p_value,med_time_high,med_time_low,fig = survival_plots(gr1_status,gr1_times,gr2_status,gr2_times,sig_gene,pval,savefile,partner)
    return(p_value,med_time_high,med_time_low,fig)

def survival_plots(event_high,time_high,event_low,time_low,gene,pval,savefile,partner):
    #former survival_analysis
    max_time = np.max(time_high+time_low)
    kmf_high = KaplanMeierFitter()
    try:
        kmf_high.fit(time_high, event_observed=event_high, label="HIGH",alpha = 1)
        med_time_high = median_survival_times(kmf_high.survival_function_)
    except ValueError:
        med_time_high = "NaN"
    kmf_low = KaplanMeierFitter()
    try:
        kmf_low.fit(time_low, event_low, label="LOW",alpha=1)
        med_time_low = median_survival_times(kmf_low.survival_function_)
    except ValueError:
        med_time_low = "NaN"

    try:
        A = len(time_high)
        B = len(time_low)
        C = len(event_high)
        D = len(event_low)
        
        if A > 0 and B > 0 and C > 0 and D >0:
            results = logrank_test(time_high, time_low, event_observed_A=event_high, event_observed_B=event_low)
            p_value = results.p_value
        else: p_value = "NA"
    except ValueError:
        p_value = "NA"
    
    fig = plt.figure(figsize=(3,3));
    plt.close(fig);
    
    if 1 > 0:#p_value < pval:
        #if round(p_value,2) <= pval and  med_time_high < med_time_low:
        if 1 > 0:
            
            plt.figure(figsize=(4, 4));
            palette={1:'b', 2:'r', 3:'g'}
            
            ax = kmf_high.plot(show_censors=False, censor_styles={'ms': 6, 'marker': 's'},c='red')
            ax = kmf_low.plot(show_censors=False, censor_styles={'ms': 6, 'marker': 's'},c='#4B73B2')#c='dodgerblue')
    
            plt.grid(False)
            
            
            ax.set_ylabel("Survival probability",fontsize=18)
            ax.set_xlabel("Years",fontsize=18)
            ax.tick_params(axis='x', labelsize=14)
            ax.tick_params(axis='y', labelsize=14)
            ax.set_yticks([0,0.2,0.4,0.6,0.8,1])
            ax.text(0,0.1,"P-value: "+str(np.round(p_value,3)),fontsize=10)
            if partner.strip()!='':
                ax.set_title(partner+":"+gene,fontsize=18)
            else:
                ax.set_title(gene,fontsize=18)
            #if savefile == 1:
            #    file = gene+"_survival.png"
            #    ax.get_figure().savefig(file,format='png',dpi=600)
            
            '''
            fig,ax = plt.subplots(figsize=(4,4))
            kmf_high.survival_function_.plot(ax=ax,color='red')
            kmf_low.survival_function_.plot(ax=ax)
            ax.set_title(gene)
            ax.set_xlabel("Years")
            ax.set_ylabel("Survival probability")
            ax.text(0,0.1,"P-value = "+str(np.round(pval,3)))
            '''
            fig = ax.get_figure();
            plt.close();
    return (p_value,med_time_high,med_time_low,fig)

def get_mut_survival_data(clinical_f,mut_samples):
    cdr_df = pd.read_csv(clinical_f,sep='\t')
    mut_survive_df = cdr_df.loc[cdr_df['bcr_patient_barcode'].isin(mut_samples)]
    mut_survive_df.reset_index(inplace=True,drop=True)
    to_drop = []
    for i in mut_survive_df.index:
        t = str(mut_survive_df.loc[i]['OS.time'])
        if t == 'nan':
            to_drop.append(i)
    mut_survive_df = mut_survive_df.drop(to_drop,axis=0)
    mut_survive_df.reset_index(inplace=True,drop=True)
    return mut_survive_df

def get_high_low_samples(df1,binding_partner,percent_cutoff):
    df1 = df1.transpose()
    cols = df1.columns.values.tolist()
    cols1 = [x.split("|")[0].upper() for x in cols]
    index = cols1.index(binding_partner)
    col_name = cols[index]
    
    df1.sort_values(by=[col_name],ascending = False,inplace = True)
    df1 = df1.transpose()
    samples = df1.columns.values.tolist()
    num_cols = len(samples)
    number_cutoff = int(percent_cutoff*num_cols/100)
    
    high_samples = samples[:number_cutoff]
    low_samples = samples[-number_cutoff:]
    
    df_high = df1[high_samples]
    df_low = df1[low_samples]
    return df_high, df_low, col_name#high_samples,low_samples
def get_ligand_information_in_clinical(clin_genes):
    genes = clin_genes #sig_genes
    progressbar = IntProgress(min=0, max=len(genes))
    display(progressbar)

    out1 = widgets.Output(layout=Layout(height='300px'))

    vbox = VBox([out1])
    display(vbox)
    ligands_df = show_ligands(out1,progressbar,genes)
    return(ligands_df)


def get_sample_info(patient_id):
    '''
    Shows a table with patient information derived from the NCI GDC Portal
    
    patient_id : str
                 TCGA patient ID, e.g. "TCGA-EE-A3AG"
    
    '''
    cases_endpt = 'https://api.gdc.cancer.gov/cases'
    filt = {
    
            "op":"=",
            "content":{
                "field":"cases.submitter_id",
                "value":[patient_id]
            }
        }
   
    params = {'filters':json.dumps(filt),"expand":"diagnoses,demographic"}
    
    response = requests.get(cases_endpt, params = params)
    out = response.json()
    dict_sample_info = {}
    flag = 0
    try:
        dict_sample_info['uuid'] = out['data']['hits'][0]['id']
        flag = 0
    except(IndexError):
        flag = 1
    if flag == 0:
        try:
            dict_sample_info['uuid'] = out['data']['hits'][0]['id']
        except(TypeError):
            dict_sample_info['uuid'] = "---"
        try:
            dict_sample_info['patient barcode'] = out['data']['hits'][0]['submitter_id']
        except(TypeError):
            dict_sample_info['patient barcode'] = "---"
        try:
            dict_sample_info['disease type'] = out['data']['hits'][0]['disease_type']
        except(TypeError):
            dict_sample_info['disease type'] = "---"
        try:
            dict_sample_info['tissue or organ of origin'] = out['data']['hits'][0]['diagnoses'][0]['tissue_or_organ_of_origin']
        except(TypeError):
            dict_sample_info['tissue or organ of origin'] = "---"
        try:
            age = out['data']['hits'][0]['diagnoses'][0]['age_at_diagnosis']*0.0027378507871321013
            years = int(age)
            days = round((age-years)*365.25)
            dict_sample_info['age at diagnosis'] = str(years)+" years "+str(days)+" days"
        except(TypeError):
            dict_sample_info['age at diagnosis'] = "---"
        try:
            dict_sample_info['ajcc pathologic stage'] = out['data']['hits'][0]['diagnoses'][0]['ajcc_pathologic_stage']
        except(TypeError):
            dict_sample_info['ajcc pathologic stage'] = "---"
        try:
            dict_sample_info['race'] = out['data']['hits'][0]['demographic']['race']
        except(TypeError):
            dict_sample_info['race'] = "---"
        try:
            dict_sample_info['gender'] = out['data']['hits'][0]['demographic']['gender']
        except(TypeError):
            dict_sample_info['gender'] = "---"
        try:
            dict_sample_info['ethnicity'] = out['data']['hits'][0]['demographic']['ethnicity']
        except(TypeError):
            dict_sample_info['ethnicity'] = "---"
        try:
            dict_sample_info['vital status'] = out['data']['hits'][0]['demographic']['vital_status']
        except(TypeError):
            dict_sample_info['vital status'] = "---"
        tbl = []
        for k in dict_sample_info.keys():
            if k != 'uuid':
                tbl.append([k,dict_sample_info[k]])
            else:
                uuid_link = '<a href="https://portal.gdc.cancer.gov/cases/'+dict_sample_info[k]+'">'+dict_sample_info[k]+'</a>'
                tbl.append([k,uuid_link])
    
        tbl = tabulate.tabulate(tbl, tablefmt='unsafehtml')
        display(tbl)    
        print("\nThe data is derived from the NCI GDC Data Portal: https://portal.gdc.cancer.gov")
    else:
        print("TCGA patient ID",patient_id,"not found.")
        

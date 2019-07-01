# Jérémie Kalfon
# for BroadInsitute
# in 2019

from __future__ import print_function

import pdb
import pandas as pd
from taigapy import TaigaClient
tc = TaigaClient()
from bokeh.palettes import viridis
import bokeh
from bokeh.resources import CDN
import numpy as np
from bokeh.plotting import *
from bokeh.models import HoverTool

def filterProteinCoding(listofgenes, idtype='ensembl_gene_id'):
  # idtype can be of "symbol","uniprot_ids","pubmed_id","ensembl_gene_id","entrez_id","name"
  tokeep = []
  b=0
  gene_mapping = tc.get(name='hgnc-87ab', file='hgnc_complete_set-2018q3')
  for i, val in enumerate(listofgenes):
      val = val.split(".")[0]
      a = gene_mapping["locus_group"][gene_mapping[idtype]==val].values
      if len(a) > 0:
          if a[0] == "protein-coding gene":
              tokeep.append(i)
      else:
          b+=1
  print(str(b))
  return(tokeep)

def convertGenes(listofgenes, from_idtype = "ensembl_gene_id", to_idtype = "symbol"):
  # idtype can be of "symbol","uniprot_ids","pubmed_id","ensembl_gene_id","entrez_id","name"
  gene_mapping = tc.get(name='hgnc-87ab', file='hgnc_complete_set-2018q3')
  not_parsed = []
  renamed = []
  b=0
  for i, val in enumerate(listofgenes):
    if from_idtype == "ensembl_gene_id":
      val = val.split(".")[0]
      a = gene_mapping[to_idtype][gene_mapping[from_idtype]==val].values
      if len(a) > 0:
        renamed.append(a[0])
      else:
          b+=1
          not_parsed.append(i)
  print(str(b) + " could not be parsed... we don't have all genes already")
  return(renamed, not_parsed)


# do deseq data on different comparison
# do volcano plot
# do a plot of variation for genes of interest to not others. 
def scatter(data, labels=None, colors=None, importance=None):
  TOOLS="hover,crosshair,pan,wheel_zoom,zoom_in,zoom_out,box_zoom,undo,redo,reset,save,box_select,lasso_select,"
  
  radi = 5
  col = viridis(len(set(colors))) if colors is not None else '#29788E' #(viridis1)
  alpha = 0.8
  radii = []
  fill_alpha = []
  cols = []
  for i in range(data.shape[0]):
    radii.append(radi if importance is None else radi/2 + importance[i]*30)
    fill_alpha.append(alpha if importance is None else alpha-(0.2*importance[i]))
    cols.append(col[0] if colors is None else col[int(colors[i])])
  source = ColumnDataSource(data=dict(
    x=data[:,0],
    y=data[:,1],
    labels=labels,
    fill_color=cols, 
    fill_alpha=fill_alpha,
    radius=radii
  ))
  TOOLTIPS = [
    ("name", "@labels"),
    ("(x,y)", "(@x, @y)"),
  ]
  p = figure(tools=TOOLS,tooltips=TOOLTIPS)
  p.circle('x','y', fill_color='fill_color', 
    fill_alpha='fill_alpha',
    radius='radius', source=source)
  show(p)

def volcano(data,genenames=None,tohighlight=None, tooltips = [('ext_gene', '@ext_gene')], 
  title="volcano plot", xlabel='log-fold change', ylabel='-log(Q)'):
  """A function to plot the bokeh single mutant comparisons."""
  # Make the hover tool
  # data should be df gene*samples + genenames

  to_plot_not, to_plot_yes = selector(dfBetaA)
  hover = bokeh.models.HoverTool(tooltips=tooltips,
                                 names=['circles'])

  # Create figure
  p = bokeh.plotting.figure(title=title, plot_width=650, 
                            plot_height=450)

  p.xgrid.grid_line_color = 'white'
  p.ygrid.grid_line_color = 'white'
  p.xaxis.axis_label = xlabel
  p.yaxis.axis_label = ylabel

  # Add the hover tool
  p.add_tools(hover)
  p = add_points(p, to_plot_not, 'b', 'qval', 'se_b', color='#1a9641')
  p = add_points(p, to_plot_yes, 'b', 'qval', 'se_b', color='#fc8d59', alpha=0.6, outline=True)
  html = file_html(p, CDN, "my plot")
  HTML(html)

  return p


def add_points(p, df1, x, y, se_x, color='blue', alpha=0.2, outline=False):
  # Define colors in a dictionary to access them with
  # the key from the pandas groupby funciton.
  df = df1.copy()
  transformed_q = -df[y].apply(np.log10)
  df['transformed_q'] = transformed_q

  source1 = bokeh.models.ColumnDataSource(df)

  # Specify data source
  p.circle(x=x, y='transformed_q', size=7,
           alpha=alpha, source=source1,
           color=color, name='circles')
  if outline:
      p.circle(x=x, y='transformed_q', size=7,
               alpha=1,
               source=source1, color='black',
               fill_color=None, name='outlines')

  # prettify
  p.background_fill_color = "#DFDFE5"
  p.background_fill_alpha = 0.5
  
  return p

def extract(df):
  """A function to separate tfs from everything else"""
  sig = (df.qval < 0.1)# & (dfBetaA.b.abs() > 0.5)
  not_tf = (~df.target_id.isin(tf.target_id))
  is_tf = (df.target_id.isin(tf.target_id))
  to_plot_not = df[sig & not_tf]
  to_plot_yes = df[sig & is_tf]
  return to_plot_not, to_plot_yes

# What pops up on hover?



###########################################################
#
# PYDESEQ
#
##################################################################

# import rpy2.robjects as robjects
# from rpy2.robjects import pandas2ri, Formula
# pandas2ri.activate()
# from rpy2.robjects.packages import importr
# deseq = importr('DESeq2')
# '''
# Adopted from: https://stackoverflow.com/questions/41821100/running-deseq2-through-rpy2
# '''

# to_dataframe = robjects.r('function(x) data.frame(x)')

# class py_DESeq2:
#   '''
#   DESeq2 object through rpy2
#   input:
#   count_matrix: should be a pandas dataframe with each column as count, and a id column for gene id
#       example:
#       id    sampleA    sampleB
#       geneA    5    1
#       geneB    4    5
#       geneC    1    2
#   design_matrix: an design matrix in the form of pandas dataframe, see DESeq2 manual, samplenames as rownames
#               treatment
#   sampleA1        A
#   sampleA2        A
#   sampleB1        B
#   sampleB2        B
#   design_formula: see DESeq2 manual, example: "~ treatment""
#   gene_column: column name of gene id columns, exmplae "id"
#   '''
#   def __init__(self, count_matrix, design_matrix, design_formula, gene_column='id'):
#     try:
#       assert gene_column in count_matrix.columns, 'Wrong gene id column name'
#       gene_id = count_matrix[gene_column]
#     except AttributeError:
#       sys.exit('Wrong Pandas dataframe?')

#     self.dds = None
#     self.deseq_result = None
#     self.resLFC = None
#     self.comparison = None
#     self.normalized_count_matrix = None
#     self.gene_column = gene_column
#     self.gene_id = count_matrix[self.gene_column]
#     self.count_matrix = pandas2ri.py2ri(count_matrix.drop(gene_column,axis=1))
#     self.design_matrix = pandas2ri.py2ri(design_matrix)
#     self.design_formula = Formula(design_formula)


#   def run_deseq(self, **kwargs):
#     self.dds = deseq.DESeqDataSetFromMatrix(countData=self.count_matrix, 
#                                     colData=self.design_matrix,
#                                     design=self.design_formula)
#     self.dds = deseq.DESeq(self.dds, **kwargs)
#     self.normalized_count_matrix = deseq.counts(self.dds, normalized=True)

#   def get_deseq_result(self, **kwargs):

#     self.comparison = deseq.resultsNames(self.dds)

#     self.deseq_result = deseq.results(self.dds, **kwargs)
#     self.deseq_result = to_dataframe(self.deseq_result)
#     self.deseq_result = pandas2ri.ri2py(self.deseq_result) ## back to pandas dataframe
#     self.deseq_result[self.gene_column] = self.gene_id.values
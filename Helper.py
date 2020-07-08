# Jeremie Kalfon
# for BroadInsitute
# in 2019

from __future__ import print_function

import pdb
import ipdb
import pandas as pd
from math import pi
import numpy as np
import itertools
import random

from taigapy import TaigaClient
tc = TaigaClient()
from IPython import get_ipython
import subprocess
from JKBio import GCPFunction as gcp
import re
import signal
import string
import sys
from PIL import Image, ImageDraw, ImageFont
import os
import json

import bokeh
from bokeh.plotting import *
from bokeh.io import output_file, show
from bokeh.transform import linear_cmap
from bokeh.util.hex import hexbin
from bokeh.models import HoverTool, CustomJS, BasicTicker, ColorBar, ColumnDataSource, LinearColorMapper, LogColorMapper, PrintfTickFormatter
from bokeh.models.widgets import TextInput
from bokeh.models.annotations import LabelSet
from bokeh.layouts import layout, widgetbox, column, row
from bokeh.resources import CDN
from bokeh.palettes import *

import matplotlib
from matplotlib import pyplot as plt
import venn as pyvenn


def fileToList(filename):
  """
  loads an input file with a\\n b\\n.. into a list [a,b,..]
  """
  with open(filename) as f:
    return [val[:-1] for val in f.readlines()]


def listToFile(l, filename):
  """
  loads a list with [a,b,..] into an input file a\\n b\\n..
  """
  with open(filename, 'w') as f:
    for item in l:
      f.write("%s\n" % item)


def dictToFile(d, filename):
  """
  turn a dict into a json file
  """
  with open(filename, 'w') as json_file:
    json.dump(d, json_file)


def fileToDict(filename):
  """
  loads a json file into a python dict
  """
  with open(filename) as f:
    data = json.load(f)
  return data


def batchMove(l, pattern=['*.', '.*'], folder='', add=''):
  """
  moves a set of files l into a folder:

  Args:
  -----
    l: file list
    pattern: if files are a set of patterns to match
    folder: folder to move file into
    add: some additional mv parameters
  """
  for val in l:
    cmd = 'mv '
    if add:
      cmd += add + ' '
    if '*.' in pattern:
      cmd += '*'
    cmd += val
    if '.*' in pattern:
      cmd += '*'
    cmd += " " + folder
    res = os.system(cmd)
    if res != 0:
      raise Exception("Leave command pressed or command failed")


def batchRename(dt, folder='', add=''):
  """
  Given a dict renames corresponding files in a folder

  Args:
  ----
    dt: dict(currentName:newName) renaming dictionnary
    folder: folder to look into
    add: some additional mv parameters
  """
  files = os.popen('ls ' + folder).read().split('\n')
  for k, val in dt.items():
    for f in files:
      if k in f:
        cmd = 'mv '
        if add:
          cmd += add + ' '
        cmd += folder
        cmd += f
        cmd += ' '
        cmd += folder
        cmd += f.replace(k, val)
        res = os.system(cmd)
        if res != 0:
          raise Exception("Leave command pressed or command failed")


def filterProteinCoding(listofgenes, from_idtype='ensembl_gene_id'):
  """
  Given a list of genes, provide the args where the genes are protein coding genes:

  This functtion will use a file in taiga, you need taigapy installed

  Args:
  -----
    listofgenes: list of genes
    from_idtype: one of "symbol","uniprot_ids","pubmed_id","ensembl_gene_id","entrez_id","name", the gene name format

  Returns:
  -------
    the args where the genes are protein coding
  """
  tokeep = []
  b = 0
  print("you need access to taiga for this (https://pypi.org/project/taigapy/)")
  gene_mapping = tc.get(name='hgnc-87ab', file='hgnc_complete_set')
  for i, val in enumerate(listofgenes):
    if from_idtype == "ensembl_gene_id":
      val = val.split(".")[0]
    elif from_idtype == "hgnc_id":
      val = "HGNC:" + str(val)
    a = gene_mapping["locus_group"][gene_mapping[from_idtype] == val].values
    if len(a) > 0:
      if a[0] == "protein-coding gene":
        tokeep.append(i)
    else:
      b += 1
  print(str(b))
  return(tokeep)


def convertGenes(listofgenes, from_idtype="ensembl_gene_id", to_idtype="symbol"):
  """
  Given a list of genes, provide the args where the genes are protein coding genes:

  This functtion will use a file in taiga, you need taigapy installed

  Args:
  -----
    listofgenes: list of genes
    from_idtype: one of "symbol","uniprot_ids","pubmed_id","ensembl_gene_id","entrez_id","name", the gene name format
    to_idtype: one of "symbol","uniprot_ids","pubmed_id","ensembl_gene_id","entrez_id","name", the gene name format

  Returns:
  -------
    1: the new names for each genes that were matched else the same name
    2: the names of genes that could not be matched
  """
  print("you need access to taiga for this (https://pypi.org/project/taigapy/)")
  gene_mapping = tc.get(name='hgnc-87ab', file='hgnc_complete_set')
  not_parsed = []
  renamed = []
  b = 0
  to = {}
  for i, val in gene_mapping.iterrows():
    to[val[from_idtype]] = val[to_idtype]
  for i, val in enumerate(listofgenes):
    if from_idtype == "ensembl_gene_id":
      val = val.split(".")[0]
    elif from_idtype == "hgnc_id":
      val = "HGNC:" + str(val)
    try:
      a = to[val]
      renamed.append(int(a) if to_idtype == 'entrez_id' else a)
    except KeyError:
      b += 1
      not_parsed.append(val)
      renamed.append(val)
  print(str(b) + " could not be parsed... we don't have all genes already")
  return(renamed, not_parsed)


def scatter(data, labels=None, title='scatter plot', showlabels=False,
            colors=None, importance=None, radi=5, alpha=0.8, **kwargs):
  """
  Makes an interactive scatter plot using Bokeh

  Args:
  -----
    data: an array-like with shape [N,2]
    labels: a list of N names for each points
    title: the plot title
    showlabels: if the labels shoul be always displayed or not (else just on hover)
    colors: a list of N integers from 0 up to 256 for the dot's colors
    importance: a list of N values to scale the size of the dots and their opacity by
    radi: the size of the dots
    alpha: the opacity of the dots
    **args: additional bokeh.figure args
  Returns:
  ------
    the bokeh object
  """
  TOOLS = "hover,crosshair,pan,wheel_zoom,zoom_in,zoom_out,box_zoom,undo,redo,reset,save,box_select,lasso_select,"

  col = viridis(len(set(colors))) if colors is not None else ['#29788E']  # (viridis1)
  radii = []
  fill_alpha = []
  cols = []
  for i in range(data.shape[0]):
    radii.append(radi if importance is None else radi / 2 + importance[i] * 30)
    fill_alpha.append(alpha if importance is None else alpha - (0.2 * importance[i]))
    cols.append(col[0] if colors is None else col[int(colors[i])])
  source = ColumnDataSource(data=dict(
      x=data[:, 0],
      y=data[:, 1],
      labels=labels if labels is not None else [''] * len(radii),
      fill_color=cols,
      fill_alpha=fill_alpha,
      radius=radii
  ))
  TOOLTIPS = [
      ("name", "@labels"),
      ("(x,y)", "(@x, @y)"),
  ]
  if showlabels:
    labels = LabelSet(x='x', y='y', text='names', level='glyph', text_font_size='9pt',
                      x_offset=5, y_offset=5, source=source, render_mode='canvas')
    p.add_layout(labels)
  p = figure(tools=TOOLS, tooltips=TOOLTIPS, title=title)
  p.circle('x', 'y', color='fill_color',
           fill_alpha='fill_alpha',
           line_width=0,
           radius = 'radius' if radi else None, source = source, kwargs)
  p.xaxis[0].axis_label=xname
  p.yaxis[0].axis_label=yname

  show(p)
  return(p)


def bigScatter(data, precomputed = False, logscale = False, features = False, colors = None, title = "BigScatter", binsize = 0.1, folder = "", showpoint = False):
  TOOLS="wheel_zoom,zoom_in,zoom_out,box_zoom,undo,redo,reset,save,box_select,lasso_select,"
  names=[("count", "@c")]
  if features:
    names.append(('features', '@features'))
  if precomputed:
    TOOLS="hover," + TOOLS
  p=figure(title = title, tools = TOOLS, tooltips = names if precomputed else None,
             match_aspect = True, background_fill_color = '#440154')
  if precomputed:
    p.hex_tile(q = "q", r = "r", size = binsize, line_color = None, source = data,
               hover_color = "pink", hover_alpha = 0.8,
               fill_color = linear_cmap('c', 'Viridis256', 0, max(data.c)) if not logscale else {'field': 'c', 'transform': LogColorMapper('Viridis256')})
  else:
    if features:
      print("we cannot yet process features on non precomputed version")
    r, bins=p.hexbin(data[:, 0], data[:, 1], line_color = None, size = binsize,
                       hover_color = "pink", hover_alpha = 0.8,
                       fill_color = linear_cmap('c', 'Viridis256', 0, None) if not logscale else {'field': 'c', 'transform': LogColorMapper('Viridis256')})
  p.grid.visible=False
  if showpoint:
    p.circle(data[:, 0], data[:, 1], color = "white", size = 1)

  if not precomputed:
    p.add_tools(HoverTool(
        tooltips=names,
        mode="mouse", point_policy="follow_mouse", renderers=[r] if not precomputed else None
    ))

  output_file(folder + title + "_hex_plot.html")

  show(p)


def CNV_Map(df, sample_order = [], title = "CN heatmaps sorted by SMAD4 loss, pointing VPS4B",
            width=900, height=400, standoff=10, y_label='', marks=[]):
  """
  create an interactive plot suited for visualizing segment level CN data for a set of samples using bokeh

  args:
  ----
    df: df['Sample' 'Start' 'End' 'Segment_Mean' 'size'] the df containing segment level copy number (can be subsetted to a specific region or genome-wide)
    sampleorder: list[Sample] <- for all samples present in the df
    title: plot title
    width: int width
    height: int height
    standoff: the space between the plot and the x axis
    y_label: the y axis label
    marks: location of lines at specific loci

  Returns:
  --------
    The bokeh object
  """
  colors = ["#75968f", "#a5bab7", "#c9d9d3", "#e2e2e2", "#dfccce", "#ddb7b1", "#cc7878", "#933b41", "#550b1d"]
  colors = RdBu[8]
  mapper = LinearColorMapper(palette=colors, low=df.Segment_Mean.min(), high=df.Segment_Mean.max())
  if len(sample_order) == 0:
    sample_order = list(set(df.Sample.tolist()))
  TOOLS = "hover,save,pan,box_zoom,reset,wheel_zoom"
  p = figure(title=title,
             y_range=(df.End.max(), df.Start.min()),
             x_range=sample_order,
             x_axis_location="above", plot_width=width, plot_height=height,
             tools=TOOLS, toolbar_location='below',
             tooltips=[('pos', '@Start, @End'), ('relative CN', '@Sample')])

  p.grid.grid_line_color = None
  p.axis.axis_line_color = None
  p.axis.major_tick_line_color = None
  p.axis.major_label_text_font_size = "5pt"
  p.axis.major_label_standoff = standoff
  p.xaxis.major_label_orientation = pi / 3
  pos = 0
  # for i,val in enumerate(historder):
  #    p.rect(x=pos,y=-7,width=len(orderbyhist[val]), height=10, fill_color=small_palettes['Viridis'][4][i])
  #    p.text(x=pos+len(orderbyhist[val])/2, y=-9, text=str(val), text_color="#96deb3")
  #    pos += len(orderbyhist[val])

  p.rect(x="Sample", y="Start", width=0.9, height="size",
         source=df.reset_index(drop=True),
         fill_color={'field': 'Segment_Mean', 'transform': mapper},
         line_color=None)

  color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="5pt",
                       ticker=BasicTicker(desired_num_ticks=len(colors)),
                       formatter=PrintfTickFormatter(format="%.2f"),
                       label_standoff=6, border_line_color=None, location=(0, 0))
  p.add_layout(color_bar, 'right')
  p.yaxis.axis_label = y_label
  # p.yaxis.major_label_overrides={20:'Centromer'}
  for val in marks:
    hline = Span(location=val, dimension='width', line_color='green', line_width=0.2)
    p.renderers.extend([hline])

  show(p)      # show the plot


def volcano(data, genenames=None, tohighlight=None, tooltips=[('gene', '@gene_id')],
            title="volcano plot", xlabel='log-fold change', ylabel='-log(Q)', maxvalue=250,
            searchbox=False, logfoldtohighlight=0.15, pvaltohighlight=0.1, showlabels=False):
  """
  Make an interactive volcano plot from Differential Expression analysis tools outputs

  Args:
  -----
    data: a df with rows genes and cols [log2FoldChange, pvalue, gene_id]
    genenames: 
    tohighlight: a list of genes to highlight in the plot
    tooltips: if user wants tot specify another bokeh tooltip
    title: plot title
    xlabel: if user wants tot specify another
    ylabel: if user wants tot specify another
    maxvalue: the max -log2(pvalue authorized usefull when managing inf vals)
    searchbox: whether or not to add a searchBox to interactively highlight genes
    minlogfold: otherwise the point is not plotted
    minpval: otherwise the point is not plotted
    logfoldtohighlight:
    pvaltohighlight:

  Returns:
  --------
    The bokeh object
  """
  to_plot_not, to_plot_yes = selector(data, tohighlight if tohighlight is not None else [], logfoldtohighlight, pvaltohighlight)
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
  p, source1 = add_points(p, to_plot_not, 'log2FoldChange', 'pvalue', color='#1a9641', maxvalue=maxvalue)
  p, source2 = add_points(p, to_plot_yes, 'log2FoldChange', 'pvalue', color='#fc8d59', alpha=0.6, outline=True, maxvalue=maxvalue)
  if showlabels:
    labels = LabelSet(x='log2FoldChange', y='transformed_q', text_font_size='7pt', text="gene_id", level="glyph",
                      x_offset=5, y_offset=5, source=source2, render_mode='canvas')
    p.add_layout(labels)
  if searchbox:
    text = TextInput(title="text", value="gene")
    text.js_on_change('value', CustomJS(
        args=dict(source=source1), code="""
      var data = source.data
      var value = cb_obj.value
      var gene_id = data.gene_id
      var a = -1
      for (i=0; i < gene_id.length; i++) {
          if ( gene_id[i]===value ) { a=i; console.log(i); data.size[i]=7; data.alpha[i]=1; data.color[i]='#fc8d59' }
      }
      source.data = data
      console.log(source)
      console.log(cb_obj)
      source.change.emit()
      console.log(source)
      """))
    p = column(text, p)
  return p


def add_points(p, df1, x, y, color='blue', alpha=0.2, outline=False, maxvalue=100):
  """parts of volcano plot"""
  # Define colors in a dictionary to access them with
  # the key from the pandas groupby funciton.
  df = df1.copy()
  transformed_q = -df[y].apply(np.log10).values
  transformed_q[transformed_q == np.inf] = maxvalue
  df['transformed_q'] = transformed_q
  df['color'] = color
  df['alpha'] = alpha
  df['size'] = 7
  source1 = bokeh.models.ColumnDataSource(df)

  # Specify data source
  p.scatter(x=x, y='transformed_q', size='size',
            alpha='alpha', source=source1,
            color='color', name='circles')
  if outline:
    p.scatter(x=x, y='transformed_q', size=7,
              alpha=1,
              source=source1, color='black',
              fill_color=None, name='outlines')

  # prettify
  p.background_fill_color = "#DFDFE5"
  p.background_fill_alpha = 0.5
  return p, source1


def selector(df, valtoextract=[], logfoldtohighlight=0.15, pvaltohighlight=0.1, minlogfold=0.15, minpval=0.1):
  """Part of Volcano plot: A function to separate tfs from everything else"""
  toshow = (df.pvalue < minpval) & (abs(df.log2FoldChange) > minlogfold)
  df = df[toshow]
  sig = (df.pvalue < pvaltohighlight) & (abs(df.log2FoldChange) > logfoldtohighlight)
  if valtoextract:
    not_tf = (~df.gene_id.isin(valtoextract))
    is_tf = (df.gene_id.isin(valtoextract))
    to_plot_not = df[~sig | not_tf]
    to_plot_yes = df[sig & is_tf]
  else:
    to_plot_not = df[~sig]
    to_plot_yes = df[sig]
  return to_plot_not, to_plot_yes

# What pops up on hover?


def plotCorrelationMatrix(data, names, colors=None, title=None, dataIsCorr=False,
                          invert=False, size=40, interactive=False, rangeto=None):
  """
  Make an interactive correlation matrix from an array using bokeh

  Args:
  -----
    data: arrayLike of int / float/ bool of size(names*val) or (names*names)
    names: list of names for each rows
    colors: list of size(names) a color for each names (good to display clusters)
    title: the plot title
    dataIsCorr: if not true, we will compute the corrcoef of the data array
    invert: whether or not to invert the matrix before running corrcoef
    size: the plot size
    interactive: whether or not to make the plot interactive (else will use matplotlib)
    rangeto: unused for now

  Returns:
    the bokeh object if interactive else None

  """
  if not dataIsCorr:
    data = np.corrcoef(np.array(data) if not invert else np.array(data).T)
  else:
    data = np.array(data)

  TOOLS = "hover,crosshair,pan,wheel_zoom,zoom_in,zoom_out,box_zoom,undo,redo,reset,save"
  xname = []
  yname = []
  color = []
  alpha = []
  if interactive:
    xname = []
    yname = []
    if rangeto is None:
      rangeto = range(len(data))
    color = []
    for i, name1 in enumerate(names):
      for j, name2 in enumerate(names):
        xname.append(name1)
        yname.append(name2)
        alpha.append(min(abs(data[i, j]), 0.9))
        if colors is not None:
          if colors[i] == colors[j]:
            color.append(Category10[10][colors[i]])
          else:
            color.append('lightgrey')
        else:
          color.append('grey' if data[i, j] > 0 else Category20[3][2])
    data = dict(
        xname=xname,
        yname=yname,
        colors=color,
        alphas=alpha,
        data=data
    )
    hover = HoverTool(tooltips=[('names: ', '@yname, @xname')])
    p = figure(title=title if title is not None else "Correlation Matrix",
               x_axis_location="above", tools=TOOLS,
               x_range=list(reversed(names)), y_range=names,
               tooltips=[('names', '@yname, @xname'), ('corr:', '@data')])

    p.plot_width = 800
    p.plot_height = 800
    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "5pt"
    p.axis.major_label_standoff = 0
    p.xaxis.major_label_orientation = np.pi / 3

    p.rect('xname', 'yname', 0.9, 0.9, source=data,
           color='colors', alpha='alphas', line_color=None,
           hover_line_color='black', hover_color='colors')
    try:
      show(p)
    except:
      show(p)
    save(p, title + '.html')

    return p  # show the plot
  else:
    plt.figure(figsize=(size, 200))
    plt.title('the correlation matrix')
    plt.imshow(data)
    plt.savefig(title + ".pdf")
    plt.show()


def venn(inp, names, title="venn"):
  """
  Plots a venn diagram using the pyvenn package

  Args:
    inp: a list of sets of values (e.g. [(1,2,3,4),(2,3),(1,3,4,5)]) 
    names: list of the name of each leaf
    title: the plot title
  """
  labels = pyvenn.get_labels(inp, fill=['number', 'logic'])
  if len(inp) == 2:
    fig, ax = pyvenn.venn2(labels, names=names)
  elif len(inp) == 3:
    fig, ax = pyvenn.venn3(labels, names=names)
  elif len(inp) == 4:
    fig, ax = pyvenn.venn4(labels, names=names)
  elif len(inp) == 5:
    fig, ax = pyvenn.venn5(labels, names=names)
  elif len(inp) == 6:
    fig, ax = pyvenn.venn6(labels, names=names)
  else:
    raise ValueError('need to be between 2 to 6')
  ax.set_title(title)
  fig.savefig(title + '.png')
  fig.show()
  plt.pause(0.1)


def grouped(iterable, n):
  """
  iterate over element of list 2 at a time python
  s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ...
  """
  it = iter(iterable)
  while True:
    chunk = tuple(itertools.islice(it, n))
    if not chunk:
      return
    yield chunk


def mergeImages(images, outputpath):
  """
  will merge a set of images in python

  Args:
  -----
    images: list of image filepath
    outputpath: where to save the resulting merger
  """
  images = list(map(Image.open, images))
  widths, heights = zip(*(i.size for i in images))

  total_width = max(widths)
  max_height = sum(heights)

  new_im = Image.new('RGB', (total_width, max_height))

  y_offset = 0
  for im in images:
    new_im.paste(im, (0, y_offset))
    y_offset += im.size[1]

  new_im.save(outputpath)


def addTextToImage(image, text, outputpath, xy=(0, 0), color=(0, 0, 0), fontSize=64):
  """
  will add some text to an image in python

  Args:
  ----
    image: the image filepath
    text: the text to write
    outputpath: the location of the resulting image
    xy: the location of the text
    color: tuple(a,b,c) a tuple of 3 ints between 0 and 256
    fontSize: an int for the font size
  """
  # adds black text to the upper left by default, Arial size 64
  img = Image.open(image)
  draw = ImageDraw.Draw(img)
  # the below file path assumes you're operating macOS
  font = ImageFont.truetype("/Library/Fonts/Arial.ttf", fontSize)
  draw.text(xy, text, color, font=font)
  img.save(outputpath)


def overlap(interval1, interval2):
  """
  Given [0, 4] and [1, 10] returns [1, 4]
  Given [0, 4] and [8, 10] returns False
  """
  if interval2[0] <= interval1[0] <= interval2[1]:
    start = interval1[0]
  elif interval1[0] <= interval2[0] <= interval1[1]:
    start = interval2[0]
  else:
    return False

  if interval2[0] <= interval1[1] <= interval2[1]:
    end = interval1[1]
  elif interval1[0] <= interval2[1] <= interval1[1]:
    end = interval2[1]
  else:
    return False

  return (start, end)


def union(interval1, interval2):
  """
  Given [0, 4] and [1, 10] returns [0, 10]
  Given [0, 4] and [8, 10] returns False
  """
  if interval1[0] <= interval2[0] <= interval1[1]:
    start = interval1[0]
    end = interval1[1] if interval2[1] <= interval1[1] else interval2[1]
  elif interval1[0] <= interval2[1] <= interval1[1]:
    start = interval2[0] if interval2[0] <= interval1[0] else interval1[0]
    end = interval1[1]
  else:
    return False
  return (start, end)


def nans(df): return df[df.isnull().any(axis=1)]


def createFoldersFor(filepath):
  """
  will recursively create folders if needed until having all the folders required to save the file in this filepath
  """
  prevval = ''
  for val in filepath.split('/')[:-1]:
    prevval += val + '/'
    if not os.path.exists(prevval):
      os.mkdir(prevval)


def randomString(stringLength=6, stype='all', withdigits=True):
  """
  Generate a random string of letters and digits

  Args:
  -----
    stringLength: the amount of char
    stype: one of lowercase, uppercase, all
    withdigits: digits allowed in the string?

  Returns:
  -------
    the string
  """
  if stype == 'lowercase':
    lettersAndDigits = ascii_lowercase
  elif stype == 'uppercase':
    lettersAndDigits = ascii_uppercase
  else:
    lettersAndDigits = string.ascii_letters
  if withdigits:
    lettersAndDigits += string.digits
  return ''.join(random.choice(lettersAndDigits) for i in range(stringLength))


def pdDo(df, op="mean", of="value1", over="value2"):
  """
  apply a function to a panda dataframe WIP
  """
  df = df.sort_values(by=over)
  index = []
  data = df.iloc[0, of]
  ret = []
  prev = df.iloc[0, over]
  j = 0
  for k, val in df.iloc[1:].iterrows():
    if val[over] == prev:
      data.append(val[of])
    else:
      if of == "mean":
        ret[j] = np.mean(data)
      elif of == "sum":
        ret[j] = np.sum(data)
      elif of == "max":
        ret[j] = np.max(data)
      elif of == "min":
        ret[j] = np.min(data)
      index.append(k)
      j += 1
      data = [val[of]]
  return index, ret


def parrun(cmds, cores, add=[]):
  """
  runs a set of commands in parallel using the "&" command

  Args:
  -----
    cmds: the list of commands
    cores: number of parallel execution
    add: an additional list(len(cmds)) of command to run in parallel at the end of each parallel run
  """
  count = 0
  exe = ''
  if len(add) != 0 and len(add) != len(cmds):
    raise ValueError("we would want them to be the same size")
  else:
    addexe = ''
  for i, cmd in enumerate(cmds):
    count += 1
    exe += cmd
    if len(add) != 0:
      addexe += add[i]
    if count < cores and i < len(cmds) - 1:
      exe += ' & '
      if len(add) != 0:
        addexe += ' & '
    else:
      count = 0
      res = subprocess.run(exe, capture_output=True, shell=True)
      if res.returncode != 0:
        raise ValueError('issue with the command: ' + str(res.stderr))
      exe = ''
      if len(add) != 0:
        res = subprocess.run(addexe, capture_output=True, shell=True)
        if res.returncode != 0:
          raise ValueError('issue with the command: ' + str(res.stderr))
        addexe = ''


def askif(quest):
  """
  asks a y/n question to the user about something and returns true or false given his answer
  """
  print(quest)
  inp = input()
  if inp in ['yes', 'y', 'Y', 'YES', 'oui', 'si']:
    return 1
  elif inp in ['n', 'no', 'nope', 'non', 'N']:
    return 0
  else:
    return askif('you need to answer by yes or no')


def inttodate(i, lim=1965, unknown='U', sep='-', order="asc", startsatyear=0):
  """
  transforms an int representing days into a date

  Args:
  ----
    i: the int
    lim: the limited year below which we have a mistake
    unknown: what to return when unknown (date is bellow the limited year)
    sep: the sep between your date (e.g. /, -, ...)
    order: if 'asc', do d,m,y else do y,m,d
    startsatyear: when is the year to start counting for this int

  Returns:
    the date or unknown
  """
  a = int(i // 365)
  if a > lim:
    a = str(a+startsatyear)
    r = i % 365
    m = str(int(r // 32)) if int(r // 32) > 0 else str(1)
    r = r % 32
    d = str(int(r)) if int(r) > 0 else str(1)
  else:
    return unknown
  return d + sep + m + sep + a+ if order == "asc" else a + sep + m + sep + d


def datetoint(dt, split='-', unknown='U', order="des"):
  """
  same as inttodate but in the opposite way;

  starts at 0y,0m,0d

  dt: the date string
  split: the splitter in the string (e.g. /,-,...)
  unknown: maybe the some dates are 'U' or 0 and the program will output 0 for unknown instead of crashing
  order: if 'asc', do d,m,y else do y,m,d

  Returns:
    the date
  """
  arr = np.array(dt[0].split(split) if dt[0] != unknown else [0, 0, 0]).astype(int)
  if len(dt) > 1:
    for val in dt[1:]:
      arr = np.vstack((arr, np.array(val.split(split) if val != unknown else [0, 0, 0]).astype(int)))
    arr = arr.T
  res = arr[2] * 365 + arr[1] * 31 + arr[0] if order == "asc" else arr[0] * 365 + arr[1] * 31 + arr[2]
  return [res] if type(res) is np.int64 else res


def getBamDate(bams, split='-', order="des", unknown='U'):
  """
  from bam files (could be in a google bucket) returns their likely sequencing date if available in the header

  Args:
  -----
    bams: the bams file|bucket paths 
    split: the splitter in the output date
    unknown: maybe the some dates can't be found the program will output unknown for them
    order: if 'asc', do d,m,y else do y,m,d

  Returns:
  -------
    a list of likely dates or [unknown]s
  """
  DTs = []
  for i, bam in enumerate(bams):
    print(i / len(bams), end='\r')
    data = os.popen('export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`\
       && samtools view -H ' + bam + ' | grep "^@RG"')
    if data == signal.SIGINT:
      print('Awakened')
      break
    else:
      res = data.read()
      dt = re.findall("(?<=\tDT:).+?\t", res)
    if len(dt) > 1:
      arr = np.array(dt[0].split('T')[0].split(split)).astype(int)
      for val in dt[1:]:
        arr = np.vstack((arr, np.array(val.split('T')[0].split(split)).astype(int)))
      arr = arr.T
      i = arr[0] * 365 + arr[1] * 31 + arr[2] if order == "asc" else arr[2] * 365 + arr[1] * 31 + arr[0]
      DTs.append(dt[np.argsort(i)[0]].split('T')[0])
    elif len(dt) == 1:
      DTs.append(dt[0].split('T')[0])
    else:
      DTs.append(unknown)
  return DTs


def getSpikeInControlScales(refgenome, fastq=None, fastQfolder='', mapper='bwa', pairedEnd=False, cores=1,
                            pathtosam='samtools', pathtotrim_galore='trim_galore', pathtobwa='bwa',
                            totrim=True, tomap=True, tofilter=True, results='res/', toremove=False):
  """
  Will extract the spikeInControls from a fastq file (usefull for, let say ChIPseq data with spike ins)

  Count based sequencing data is not absolute and will be normalized as each sample will be sequenced at a specific depth. 
  To figure out what was the actual sample concentration, we use Spike In control
  You should have FastQfolder/[NAME].fastq & BigWigFolder/[NAME].bw with NAME being the same for the same samples


  Args:
  -----
  refgenome: str the file path to the indexed reference genome
  FastQfolder: str the folder path where the fastq files are stored (should be named the same as files in BigWigFolder)
  BigWigFolder: str the folder path where the bigwig files are stored (should be named the same as files in FastQfolder)
  mapper: str flag to 'bwa', ...
  pairedEnd: Bool flat to true for paired end sequences. if true, You should have FastQfolder/[NAME]_1|2.fastq

  Returns:
  --------
  dict(file,float) the scaling factor dict

  """
  if len(fastQfolder) > 0:
    print('using all files from folder')
    fastqs = os.listdir(fastQfolder)
    fastqs = [i for i in fastqs if '.fq.gz' == i[-6:] or '.fastq.gz' == i[-9:]]
    fastqs.sort()
    if pairedEnd and (tomap or totrim):
      print("need to be name_*1, name_*2")
      fastqs = [i for i in grouped(fastqs, 2)]
  elif fastq is None:
    raise Error('you need input files')
  else:
    if type(fastq) is list:
      print('your files need to be all in the same folder')
      fastQfolder = '/'.join(fastq[0].split('/')[:-1]) + '/'
      if not totrim and not tomap:
        fastqs = [f.split('/')[-1] for f in fastq]
      else:
        print("need to be name_*1, name_*2")
        fastqs = [[f[0].split('/')[-1], f[1].split('/')[-1]] for f in grouped(fastq, 2)]
    else:
      fastQfolder = '/'.join(fastq.split('/')[:-1]) + '/'
      fastqs = [fastq.split('/')[-1]]
  print(fastqs)
  if not totrim:
    print("you need to have your files in the " + results + " folder")
  if totrim and tomap:
    print("\n\ntrimming\n\n")
    if pairedEnd:
      cmds = []
      rm = []
      for file in fastqs:
        cmd = pathtotrim_galore + ' --paired --fastqc --gzip ' + fastQfolder + file[0] + ' ' + fastQfolder + file[1] + " -o " + results
        if toremove:
          rm.append('rm ' + fastQfolder + file[0] + ' ' + fastQfolder + file[1])
        cmds.append(cmd)
      print(cmds)
      parrun(cmds, cores, add=rm)
      fastqs = [[file[0].split('.')[0] + '_val_1.fq.gz', file[1].split('.')[0] + '_val_2.fq.gz'] for file in fastqs]
  if tomap:
    print("\n\nmapping\n\n")
    if pairedEnd:
      cmds = []
      rm = []
      for file in fastqs:
        cmd = pathtobwa + ' mem ' + refgenome + ' ' + results + file[0] + ' ' + results +\
            file[1] + ' | ' + pathtosam + ' sort - -o ' + results + file[0].split('.')[0] + '.sorted.bam'
        if toremove:
          rm.append('rm ' + results + file[0] + ' ' + results + file[1])
        cmds.append(cmd)
      parrun(cmds, cores, add=rm)
      fastqs = [file[0].split('.')[0] + '.sorted.bam' for file in fastqs]

  if tofilter:
    print("\n\nfiltering\n\n")
    cmds = []
    rm = []
    parrun([pathtosam + ' index ' + results + file.split('.')[0] + '.sorted.bam' for file in fastqs], cores)
    parrun([pathtosam + ' flagstat ' + results + file.split('.')[0] + '.sorted.bam > ' + results + file.split('.')[0] + '.sorted.bam.flagstat' for file in fastqs], cores)
    parrun([pathtosam + ' idxstats ' + results + file.split('.')[0] + '.sorted.bam > ' + results + file.split('.')[0] + '.sorted.bam.idxstat' for file in fastqs], cores)
    fastqs = [file.split('.')[0] + '.sorted.bam' for file in fastqs]
  else:
    print("files need to be named: NAME.sorted.bam")
    fastqs = [file for file in fastqs if '.sorted.bam' == file[-11:]]
  mapped = {}
  norm = {}
  unique_mapped = {}
  print("\n\ncounting\n\n")
  for file in fastqs:
    mapped[file.split('.')[0]] = int(os.popen(pathtosam + ' view -c -F 0x004 -F 0x0008 -f 0x001 -F 0x0400 -q 1 ' + results +
                                              file + ' -@ ' + str(cores)).read().split('\n')[0])
   # unique_mapped[file.split('.')[0]] = int(re.findall("Mapped reads: (\d+)", os.popen('bamtools stats -in '+results +
    #                                             file + '.sorted.bam').read())[0])
  nbmapped = np.array([i for i in mapped.values()])
  nbmapped = np.sort(nbmapped)[0] / nbmapped.astype(float)
  for i, val in enumerate(mapped.keys()):
    norm[val] = nbmapped[i]
  return norm, mapped,  # unique_mapped


def changeToBucket(samples, gsfolderto, values=['bam', 'bai'], catchdup=False):
  """
  """
  # to do the download to the new dataspace
  for i, val in samples.iterrows():
    ran = randomString(6, 'underscore', withdigits=False)
    for ntype in values:
      name = val[ntype].split('/')[-1] if catchdup else ran + '_' + val[ntype].split('/')[-1]
      if not gcp.exists(gsfolderto + val[ntype].split('/')[-1]) or not catchdup:
        cmd = 'gsutil cp ' + val[ntype] + ' ' + gsfolderto + name
        res = subprocess.run(cmd, shell=True, capture_output=True)
        if res.returncode != 0:
          raise ValueError(str(res.stderr))
      else:
        print(name + ' already exists in the folder: ' + gsfolderto)
        print(gcp.lsFiles([gsfolderto + name], '-la'))
      samples.loc[i, ntype] = gsfolderto + name
  return samples


def GSEAonExperiments(data, experiments, res={}, savename='', scaling=[], geneset='GO_Biological_Process_2015',
                      cores=8, cleanfunc=lambda i: i.split('(GO')[0]):
  """

  Will run GSEA on a set of experiment
  
  Args:
  -----
    data: a pandas.df rows: gene counts; columns: [experimentA_X,..., experimentD_X..., control_X] where X is the replicate number
    experiments: a list of experiment names (here experimentA,.. experimentD)
    scaling: a dict(experiment:(mean,std)) of scaling factors and their associated standard error for each experiments
    res: you can provide a dict containing results from
    savename: if you want to save the plots as pdfs, provides a location/name
    geneset: the geneset to run it on. (can be a filepath to your own geneset)
    cores: to run GSEA on
    cleanfunc: a func applied to the names of the gene sets to change it in some way (often to make it more readable)
  Returns
  -------
    plots the results
    1: returns a matrix with the enrichment for each term for each experiment
    2: returns a dict(experiment:pd.df) with dataframe being the output of GSEA (with pvalues etc..) for each experiments
  """
  for i, val in enumerate(experiments):
    print(val)
    totest = data[[v for v in data.columns[:-1] if val + '-' in v or 'AAVS1' in v]]
    cls = ['Condition' if val + '-' in v else 'DMSO' for v in totest.columns]
    if scaling:
      if abs(scaling[val.split('_')[1]][0]) > scaling[val.split('_')[1]][1]:
        print("rescaling this one")
        cols = [i for i in totest.columns if val + '-' in i]
        totest[cols] = totest[cols] * (2**scaling[val.split('_')[1]][0])
    if val in res:
      print(val + " is already in set")
      continue
    res[val] = gseapy.gsea(data=totest, gene_sets=geneset,
                           cls=cls, no_plot=False, processes=cores)
    res[val].res2d['Term'] = [i for i in res[val].res2d.index]
    for i, v in res.items():
      res[i].res2d['Term'] = [cleanfunc(i) for i in v.res2d['Term']]
    plt.figure(i)
    sns.barplot(data=res[val].res2d.iloc[:25], x="es", y="Term",
                hue_order="geneset_size").set_title(val)
  a = set()
  for k, val in res.items():
    a.update(set(val.res2d.Term))
  a = {i: [0] * len(res) for i in a}
  for n, (k, val) in enumerate(res.items()):
    for i, v in val.res2d.iterrows():
      a[v.Term][n] = v.es
  pres = pd.DataFrame(a, index=res.keys())
  a = sns.clustermap(figsize=(25, 20), data=res, vmin=-1, vmax=1, yticklabels=res.index, cmap=plt.cm.RdYlBu)
  b = sns.clustermap(-res.T.corr(), cmap=plt.cm.RdYlBu, vmin=-1, vmax=1)
  if savename:
    res.to_csv(savename + ".csv")
    a.savefig(savename + "_genesets.pdf")
    b.savefig(savename + "_correlation.pdf")
  return pres, res


def runERCC(ERCC, experiments, featurename="Feature", issingle=False, dilution=1 / 100,
            name="RNPv2", spikevol=1, control="AAAVS1", fdr=0.1, totalrnamass=0.5):
  """
  Runs the ERCC dashboard Rpackage on your notebook

  you will need to run this function from ipython and to have the R package erccdashboard installed
  
  Args:
  ----
    ERCC: a pandas.df rows: ERCC counts columns: [experimentA_X,..., experimentD_X..., control_X] where X is the replicate number
    experiments: a list of experiment names (here experimentA,.. experimentD)
    featurename: columns where the ERCC pseudo gene names are stored
    issingle: ERCC parameters to choose between Single and RatioPair
    dilution: ERCC dilution parameter
    name: the name of the experiment set
    spikevol: ERCC spikevol parameter
    control: the control name (here control)
    fdr: ERCC fdr parameter
    totalrnamass: ERCC totalrnamass parameter
  
  Returns:
  -------
    a dict(experimentName:(val, ste)) a dict containing the scaling factor and its standard error for each experiment

  Raises:
  ------
    RuntimeError: if you are not on ipython
  """
  try:
    ipython = get_ipython()
  except:
    raise RuntimeError('you need to be on ipython')
  ipython.magic("load_ext rpy2.ipython")
  ipython.magic("R library('erccdashboard')")
  ipython.magic("R datType = 'count'")  # "count" for RNA-Seq data, "array" for microarray data
  ipython.magic("R isNorm = F")  # flag to indicate if input expression measures are already normalized, default is FALSE
  ipython.magic("R -i name filenameRoot = name")  # user defined filename prefix for results files
  ipython.magic("R -i control sample2Name = control")  # name for sample 2 in the experiment
  ipython.magic("R -i issingle erccmix < - if(issingle) 'Single' else 'RatioPair'")  # name of ERCC mixture design, "RatioPair" is default
  ipython.magic("R -i dilution erccdilution = dilution")  # dilution factor used for Ambion spike-in mixtures
  ipython.magic("R -i spikevol spikeVol = spikevol")  # volume (in microliters) of diluted spike-in mixture added to total RNA mass
  ipython.magic("R -i fdr choseFDR = fdr")  # user defined false discovery rate (FDR), default is 0.05
  ipython.magic("R exDat = ''")
  ipython.magic("R - i totalrnamass totalRNAmass < - totalrnamass")

  cols = list(ERCC.columns)
  cols.sort()
  res = {}
  for val in experiments:
    d = {}
    e = 0
    c = 0
    d.update({
        featurename: 'Feature'
    })
    for i in cols[:-1]:
      if val + '-' in i:
        e += 1
        d.update({i: val.split('_')[-1] + '_' + str(e)})
      if control + "-" in i:
        c += 1
        d.update({i: control + "_" + str(c)})
    a = ERCC[list(d.keys())].rename(columns=d)
    a.to_csv('/tmp/ERCC_estimation.csv', index=None)
    val = val.split('_')[-1]

    ipython.magic("R -i val print(val)")

    ipython.magic("R print(sample2Name)")
    ipython.magic("R a <- read.csv('/tmp/ERCC_estimation.csv')")
    ipython.magic("R print(head(a))")

    try:

      ipython.magic("R - i val exDat = initDat(datType=datType, isNorm=isNorm, exTable=a,\
                                 filenameRoot=filenameRoot, sample1Name=val,\
                                 sample2Name=sample2Name, erccmix=erccmix,\
                                 erccdilution=erccdilution, spikeVol=spikeVol,\
                                 totalRNAmass=totalRNAmass, choseFDR=choseFDR)")
      ipython.magic("R exDat = est_r_m(exDat)")
      ipython.magic("R exDat = dynRangePlot(exDat)")

    except Warning:
      print("failed for " + val)
      continue
    except:
      print('worked for ' + val)

    ipython.magic("R print(summary(exDat))")
    ipython.magic("R grid.arrange(exDat$Figures$dynRangePlot)")
    ipython.magic("R grid.arrange(exDat$Figures$r_mPlot)")
    ipython.magic("R grid.arrange(exDat$Figures$rangeResidPlot)")

    ipython.magic("R -o rm rm <- exDat$Results$r_m.res$r_m.mn")
    ipython.magic("R -o se se <- exDat$Results$r_m.res$r_m.mnse")

    res[val] = (rm[0], se[0])
  for i, v in res.items():
    if abs(v[0]) > v[1]:
      print(i, v[0])
  return res


def fromGTF2BED(gtfname, bedname, gtftype='geneAnnot'):
  """
  transforms a  gtf file into a bed file

  Args:
  ----
    gtfname: filepath to gtf file
    bedname: filepath to beddfile
    gtftype: only geneAnnot for now

  Returns:
  --------
    newbed: the bedfile as a pandas.df

  """
  if gtftype == 'geneAnnot':
    gtf = pd.read_csv(gtfname, sep='\t', header=0, names=["chr", "val", "type", "start", 'stop', 'dot', 'strand', 'loc', 'name'])
    gtf['name'] = [i.split('gene_id "')[-1].split('"; trans')[0] for i in gtf['name']]
    prevname = ''
    newbed = {'chr': [], 'start': [], 'end': [], 'gene': []}
    for i, val in gtf.iterrows():
      showcount(i, len(gtf))
      if val['name'] == prevname:
        newbed['end'][-1] = val['stop']
      else:
        newbed['chr'].append(val['chr'])
        newbed['start'].append(val['start'])
        newbed['end'].append(val['stop'])
        newbed['gene'].append(val['name'])
      prevname = val['name']
    newbed = pd.DataFrame(newbed)
    newbed = newbed[~newbed.chr.str.contains('_fix')]
    newbed.to_csv(bedname + ".bed", sep='\t', index=None)
    newbed.to_csv(bedname + "_genes.bed", sep='\t', index=None)
    return newbed


def showcount(i, size):
  """
  pretty print of i/size%, to put in a for loop
  """
  print(str(1 + int(100 * (i / size))) + '%', end='\r')


def combin(n, k):
  """
  Nombre de combinaisons de n objets pris k a k
  Number of comabination of n object taken k at a time
  """
  if k > n // 2:
    k = n - k
  x = 1
  y = 1
  i = n - k + 1
  while i <= n:
    x = (x * i) // y
    y += 1
    i += 1
  return x


def mergeSplicingVariants(df, defined='.'):
  df = df.T.sort_index()
  foundpoint = False
  for i, v in enumerate(df.index.tolist()):
    if foundpoint:
      if defined in v:
        tomerge.append(v)
      else:
        if foundpoint not in df.index:
          df.loc[foundpoint] = df.loc[tomerge].sum(1)
        df = df.drop(index=tomerge)
        foundpoint = False
    elif defined in v:
      foundpoint = v.split(defined)[0]
      tomerge = [v]
  return df
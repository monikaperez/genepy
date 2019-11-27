# Jeremie Kalfon
# for BroadInsitute
# in 2019

from __future__ import print_function

import pdb
import pandas as pd
from taigapy import TaigaClient
tc = TaigaClient()
from bokeh.palettes import *
import bokeh
from bokeh.resources import CDN
import numpy as np
from bokeh.plotting import *
from bokeh.models import HoverTool

import matplotlib
import venn as pyvenn
import sys
from PIL import Image, ImageDraw, ImageFont
import os


def fileToList(filename):
  with open(filename) as f:
    return [val[:-1] for val in f.readlines()]


def filterProteinCoding(listofgenes, idtype='ensembl_gene_id'):
  # idtype can be of "symbol","uniprot_ids","pubmed_id","ensembl_gene_id","entrez_id","name"
  tokeep = []
  b = 0
  print("you need access to taiga for this (https://pypi.org/project/taigapy/)")
  gene_mapping = tc.get(name='hgnc-87ab', file='hgnc_complete_set')
  for i, val in enumerate(listofgenes):
    if idtype == "ensembl_gene_id":
      val = val.split(".")[0]
    elif idtype == "hgnc_id":
      val = "HGNC:" + str(val)
    a = gene_mapping["locus_group"][gene_mapping[idtype] == val].values
    if len(a) > 0:
      if a[0] == "protein-coding gene":
        tokeep.append(i)
    else:
      b += 1
  print(str(b))
  return(tokeep)


def convertGenes(listofgenes, from_idtype="ensembl_gene_id", to_idtype="symbol"):
  # idtype can be of "symbol","uniprot_ids","pubmed_id","ensembl_gene_id","entrez_id","name"
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
      renamed.append(int(a) if from_idtype == 'entrez_id' else a)
    except KeyError:
      b += 1
      not_parsed.append(val)
      renamed.append(val)
  print(str(b) + " could not be parsed... we don't have all genes already")
  return(renamed, not_parsed)


# do deseq data on different comparison
# do volcano plot
def scatter(data, labels=None, colors=None, importance=None, radi=5, alpha=0.8, **kargs):
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
  p = figure(tools=TOOLS, tooltips=TOOLTIPS)
  p.circle('x', 'y', color='fill_color',
           fill_alpha='fill_alpha',
           line_width=0,
           radius='radius', source=source)
  return(p)


def bar():
  data['Year'] = data['Year'].astype(str)
  data = data.set_index('Year')
  data.drop('Annual', axis=1, inplace=True)
  data.columns.name = 'Month'

  years = list(data.index)
  months = list(data.columns)

  # reshape to 1D array or rates with a month and year for each row.
  df = pd.DataFrame(data.stack(), columns=['rate']).reset_index()

  # this is the colormap from the original NYTimes plot
  colors = ["#75968f", "#a5bab7", "#c9d9d3", "#e2e2e2", "#dfccce", "#ddb7b1", "#cc7878", "#933b41", "#550b1d"]
  mapper = LinearColorMapper(palette=colors, low=df.rate.min(), high=df.rate.max())

  TOOLS = "hover,save,pan,box_zoom,reset,wheel_zoom"

  p = figure(title="US Unemployment ({0} - {1})".format(years[0], years[-1]),
             x_range=years, y_range=list(reversed(months)),
             x_axis_location="above", plot_width=900, plot_height=400,
             tools=TOOLS, toolbar_location='below',
             tooltips=[('date', '@Month @Year'), ('rate', '@rate%')])

  p.grid.grid_line_color = None
  p.axis.axis_line_color = None
  p.axis.major_tick_line_color = None
  p.axis.major_label_text_font_size = "5pt"
  p.axis.major_label_standoff = 0
  p.xaxis.major_label_orientation = pi / 3

  p.rect(x="Year", y="Month", width=1, height=1,
         source=df,
         fill_color={'field': 'rate', 'transform': mapper},
         line_color=None)

  color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="5pt",
                       ticker=BasicTicker(desired_num_ticks=len(colors)),
                       formatter=PrintfTickFormatter(format="%d%%"),
                       label_standoff=6, border_line_color=None, location=(0, 0))
  p.add_layout(color_bar, 'right')

  show(p)      # show the plot


def CNV_Map(df, sample_order=[], title="CN heatmaps sorted by SMAD4 loss, pointing VPS4B",
            width=900, height=400, standoff=10, ylabel='', marks=[]):
  """
  args:
  ----
    df: df[Sample Start End Segment_Mean size]
    sampleorder: list[Sample] <- for all samples present in the df
  """
  colors = ["#75968f", "#a5bab7", "#c9d9d3", "#e2e2e2", "#dfccce", "#ddb7b1", "#cc7878", "#933b41", "#550b1d"]
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
            title="volcano plot", xlabel='log-fold change', ylabel='-log(Q)', maxvalue=250):
  """A function to plot the bokeh single mutant comparisons."""
  # Make the hover tool
  # data should be df gene*samples + genenames
  to_plot_not, to_plot_yes = selector(data, tohighlight if tohighlight is not None else [])
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
  p = add_points(p, to_plot_not, 'log2FoldChange', 'pvalue', 'se_b', color='#1a9641', maxvalue=maxvalue)
  p = add_points(p, to_plot_yes, 'log2FoldChange', 'pvalue', 'se_b', color='#fc8d59', alpha=0.6, outline=True, maxvalue=maxvalue)
  return p


def add_points(p, df1, x, y, se_x, color='blue', alpha=0.2, outline=False, maxvalue=100):
  # Define colors in a dictionary to access them with
  # the key from the pandas groupby funciton.
  df = df1.copy()
  transformed_q = -df[y].apply(np.log10).values
  transformed_q[transformed_q == np.inf] = maxvalue
  df['transformed_q'] = transformed_q
  source1 = bokeh.models.ColumnDataSource(df)
  source2 = bokeh.models.ColumnDataSource(df)

  # Specify data source
  p.circle(x=x, y='transformed_q', size=7,
           alpha=alpha, source=source1,
           color=color, name='circles')
  if outline:
    p.circle(x=x, y='transformed_q', size=7,
             alpha=1,
             source=source2, color='black',
             fill_color=None, name='outlines')

  # prettify
  p.background_fill_color = "#DFDFE5"
  p.background_fill_alpha = 0.5

  return p


def selector(df, valtoextract):
  """A function to separate tfs from everything else"""
  sig = (df.pvalue < 0.1) & (abs(df.log2FoldChange) > 0.15)
  not_tf = (~df.gene_id.isin(valtoextract))
  is_tf = (df.gene_id.isin(valtoextract))
  to_plot_not = df[sig & not_tf]
  to_plot_yes = df[sig & is_tf]
  return to_plot_not, to_plot_yes

# What pops up on hover?


def plotCorrelationMatrix(data, names, colors=None, title=None, dataIsCorr=False,
                          invert=False, size=40, interactive=False, rangeto=None):
  """
  data arrayLike of int / float/ bool of size(names*val)
  names list like string
  colors, list like size(names)

  """
  if not dataIsCorr:
    data = np.corrcoef(np.array(data))
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

    save(p, 'temp/corrmat.html')
    return p  # show the plot
  else:
    plt.figure(figsize=(size, 200))
    plt.title('the correlation matrix')
    plt.imshow(data.T if invert else data)
    plt.savefig("temp/corrmat.pdf")
    plt.show()


def venn(inp, names):
  labels = pyvenn.get_labels(inp, fill=['number', 'logic'])
  if len(inp) == 2:
    fig, ax = pyvenn.venn2(labels, names=names)
  if len(inp) == 3:
    fig, ax = pyvenn.venn3(labels, names=names)
  if len(inp) == 4:
    fig, ax = pyvenn.venn4(labels, names=names)
  if len(inp) == 5:
    fig, ax = pyvenn.venn5(labels, names=names)
  if len(inp) == 6:
    fig, ax = pyvenn.venn6(labels, names=names)
  fig.show()


def grouped(iterable, n):
  """
  iterate over element of list 2 at a time python
  s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ...
  """
  return zip(*[iter(iterable)] * n)


def mergeImages(images, outputpath):
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


def addTextToImage(imagedir, text, outputpath, xy=(0, 0), color=(0, 0, 0), fontSize=64):
    # adds black text to the upper left by default, Arial size 64
  img = Image.open(imagedir)
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
  prevval = ''
  for val in filepath.split('/')[:-1]:
    prevval += val + '/'
    if not os.path.exists(val):
      os.mkdir(val)


def pdDo(df, op="mean", of="value1", over="value2"):
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

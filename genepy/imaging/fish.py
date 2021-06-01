import numpy as np
from scipy.spatial import distance_matrix
from genepy.utils import helper as h
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


# make a plot of averaged binned signal strength by distance from locis
def computeDistsFromClass(dots, seconddots, conds=['DMSO', 'VHL'], groupcol="group",
                          sclass='green', signal="mean_green", area="area"):
    """
    """
    dists= {}
    twodists = {}
    for val in set(dots.exp):
        for e in conds:
            d = dots[(dots.exp==val)&(dots.treat==e)]
            dist = []
            weight = []
            newdist = []
            m = seconddots[(seconddots.exp==val)&(seconddots.treat==e)]
            print(val, e)
            for i,(k, v) in enumerate(m.iterrows()):
                h.showcount(i, len(m))
                dist.append(
                    distance_matrix(d[(d['class']==sclass)&
                                      (d[groupcol]==v[groupcol])][['x', "y", "z"]].values,
                                    np.array([v[['x_mean', "y_mean", "z_mean"]]])).T[0].astype(float))
                weight.append(d[(d['class'] == sclass)&(d[groupcol]==v[groupcol])][signal])
                a = d[(d['class']==sclass) &
                      (d[groupcol]==v[groupcol])][['x',"y","z",signal,area,"m_id"]].values
                a[:,:3] = a[:,:3] - v[['x_mean', "y_mean", "z_mean"]].values
                newdist.append(a)
            twodists[val+e] = pd.DataFrame(data=np.vstack(newdist),
                                           columns=['x', 'y', 'z', signal, area, "m_id"])
            dists[val+e] = [np.hstack(dist), np.hstack(weight)]
    return  twodists, dists


def drawDots(dists, scenter=False, size=1000, zsize=1000,
             folder="", sscale=1, signal="signal",
             area="area", ** kwargs):
    """
    """
    for i, (k,a) in enumerate(dists.items()):
        a = a.copy()
        a[area] = ((a[area]/(3.14))**(1/2)).astype(float)
        if sscale:
            a[signal] /= sscale[i]
        a = a[(abs(a.x)<size*1.4) & (abs(a.y)<size*1.4) & (abs(a.z)<zsize*1.4)]

        ax = sns.scatterplot(data=a, x='x', y='y',
                            hue=signal, size=area, alpha=sscale[i] if sscale else None, **kwargs)
        plt.title(k)
        if scenter:
            ax.plot([0], [0], 'o', ms=scenter, markerfacecolor="None",
            markeredgecolor='red', markeredgewidth=1)

        ax.set_xlim((-size,size))
        ax.set_ylim((-size,size))

        plt.show()
        plt.savefig(folder+''.join(k)+"_scatter_representation_size_to_center.pdf")


def colocalize(dots, groupcol='group', zcol='z', xcol="x", ycol="y",
              default_maxdist=None, distance_scale=1.2, areacol="area",
              mergedidcol='m_id'):
    """
    """
    idcount = 0
    merged = dots.copy()
    for group in set(dots[groupcol]):
        print(group)
        # merging cells
        gdot = dots[(dots[groupcol] == group)]
        maxdist = default_maxdist if default_maxdist else np.sqrt(
            gdot[areacol].mean() / 3.14) * distance_scale

        gdot[mergedidcol]=None

        pos = gdot[[xcol, ycol, zcol]].values
        dist = distance_matrix(pos, pos)

        # closest needs to not be too far away otherwise the dot
        # is considered as finished
        for val in np.tril(dist < maxdist):

            con = np.argwhere(val > 0).T[0].tolist()
            # we get all its connections
            con_val = gdot.iloc[con]
            ids = list(set(con_val[mergedidcol]) - set([None]))
            # if connections are already connected we use this id
            if len(ids) > 0:
                def_id = ids[0]
                # for each connection, if have another id,
                # replace this id with the current one
                for i in ids:
                    con.extend(np.argwhere(
                        (gdot[mergedidcol] == i).values).T[0].tolist())
                con = list(set(con))
            # if none we create a new id
            else:
                idcount+=1
                def_id = "id_"+str(idcount)
            gdot.loc[gdot.iloc[con].index.tolist(), mergedidcol] = def_id
        #except:
        #    pdb.set_trace()
        merged.loc[gdot.index.tolist(), mergedidcol] = gdot[mergedidcol].tolist()
    return merged


def mergeAnnotated(annot, minzstack=2, groupdefault={}, todrop=[], coltocount="image",
                    id="m_id", colocName="cobinding"):
    """
    """
    annot = annot.drop(columns=todrop)
    grouping = {i: "mean" for i in annot.columns}
    grouping.update(groupdefault)
    grouping.pop(id)
    # merge into a same sample
    groups = annot.groupby(id)
    counts = groups[coltocount].count()
    merged = groups.agg(grouping)
    merged['counts'] = counts
    merged = merged[merged['counts'] >= minzstack]
    merged.columns = [i[0] if "first" in i[1]
                    else '_'.join(i) for i in merged.columns]
    #rename colors
    merged['class'] = [i[0] if len(
        i) == 1 else colocName for i in merged["class_unique"]]
    return merged.drop(columns="class_unique")

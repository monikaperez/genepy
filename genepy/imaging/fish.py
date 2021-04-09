import numpy as np
from scipy.spatial import distance_matrix
import pdb
from genepy.utils import helper as h








#####################################################
############# THESE FUNCTIONS ARE IN THE ############
############# CASE WE WORK WITH Z STACK  ############
############# THAT ARE NOT ANNOTATED TOG ############
############# ETHER                      ############
#####################################################
def _mergeCells(cells, groupcol='group', zcol='z', xcol="x", ycol="y", idcol='id'):
    idcount = 2
    cellmatch = {}
    try:
        for group in set(cells[groupcol]):
            print(group)
            # merging cells
            cellmatch[group] = {}
            prevzpos = np.array([])
            previd = []
            gcell = cells[(cells[groupcol] == group)]
            for z in set(gcell[zcol]):
                zcell = gcell[gcell[zcol] == z]
                zpos = zcell[[xcol, ycol]].values
                # continuing cells
                if prevzpos.any():
                    # computing closest centers in next zindex
                    dist = distance_matrix(prevzpos, zpos)

                    # closest needs to not be too far away otherwise the cell
                    # is considered as finished
                    argm = np.argmin(dist, 1)
                    selected = dist[range(len(prevzpos)), argm]
                    valid = selected < 500
                    ok_argm = argm[valid]
                    # if the same cell is selected and valid for two different cells:
                    for val in h.dups(ok_argm):
                        print("found DUP!")
                        pdb.set_trace()
                    for i, v in enumerate(previd):
                        if valid[i]:
                            closest_id = zcell.iloc[argm[i]][idcol]
                            cellmatch[group][v].append((z, closest_id))
                    # we set the new prev values
                    prevzpos = zpos[ok_argm]
                    previd = list(np.array(previd)[valid])
                    # the start of new cells
                    unselected = list(set(range(dist.shape[1])) - set(ok_argm))
                    if len(unselected) > 0:
                        cellmatch[group].update(
                            {i+idcount: [(z, v)] for i, v in enumerate(zcell.iloc[unselected][idcol].values)})
                        prevzpos = np.vstack([prevzpos, zpos[unselected]])
                        previd.extend(range(idcount, idcount+len(unselected)))
                        idcount += len(zpos)
                else:
                    # we create the starting cells
                    cellmatch[group] = {
                        i+idcount: [(z, v)] for i, v in enumerate(zcell[idcol].values)}
                    previd = list(range(idcount, idcount+len(zcell)))
                    idcount += len(zcell)
                    prevzpos = zpos
    except:
        pdb.set_trace()
        raise
    return cellmatch


def _relabelCells(cm, cells, dots=None, idcol='id', pidcol="parent_id", groupcol='group', zcol='z', filter_cell_without_dots=True):
    # change all dots' parent_id in this cell
    pdb.set_trace()
    for k, v in cm.items():
        gcell = cells[cells[groupcol] == k]
        gdots = dots[dots[groupcol] == k]
        for l, w in v.items():
            ndots = 0
            cids = []
            for (z, val) in w:
                # setting cell ids
                cid = gcell[(gcell[idcol] == val) &
                            (gcell[zcol] == z)].index[0]
                cells.loc[cid, idcol] = l
                cids.append(cid)

                if dots is not None:
                    #matching dot's parent ids
                    dids = gdots[(gdots[pidcol] == val) & (
                        gdots[zcol] == z)].index.tolist()
                    dots.loc[dids, pidcol] = l
                    ndots += len(dids)
            if dots is not None:
                if ndots == 0 and filter_cell_without_dots:
                    cells = cells.drop(index=cids)
                cm[k][l] = ndots
    return cm, cells, dots


def mergeCells(cells, dots=None, groupcol='group', pidcol="parent_id", zcol='z', xcol="x", ycol="y", idcol='id', filter_cell_without_dots=True):
    res = _mergeCells(cells, groupcol, zcol, xcol, ycol, idcol)
    return _relabelCells(res, cells, dots, idcol,
                         pidcol, groupcol, zcol, filter_cell_without_dots)

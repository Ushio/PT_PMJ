import hou
import json
import os
import time

# node: hou.ObjNode
def exportObjGeometry(node):
    beg_time = time.time()

    if node.renderNode() is None:
        print("    This dosen't include any renderNode.")
        return

    # print("{}".format(node.worldTransform()))
    rNode = node.renderNode() # hou.SopNode
    rGeom = rNode.geometry() # hou.Geometry

    # print("{}".format(type(rNode)))

    primitiveType = ''
    for prim in rGeom.prims():
        if prim.type() == hou.primType.Polygon:
            primitiveType = 'Polygon'
            break

    data = {
        'type' : primitiveType,
        'xform' : node.worldTransform().asTuple(),
    }

    # 'triangles'
    points = []
    for point in rGeom.points():
        p = point.position()
        points.append(p.x())
        points.append(p.y())
        points.append(p.z())

    Points = {
        'P' : points
    }
    data['Points'] = Points

    skippedPrimitive = 0
    
    PointNum = [] # for Vertices
    for prim in rGeom.prims():
        if prim.type() != hou.primType.Polygon or prim.numVertices() != 3:
            hasNonTriangles 
            continue
        
        for vertex in prim.vertices():
            PointNum.append(vertex.point().number())

    Vertices = {
        'Point Num' : PointNum
    }
    data['Vertices'] = Vertices

    if skippedPrimitive != 0:
        print("    {} prims were skipped by some reasons.".format(skippedPrimitive))

    # output Path
    filePath = os.path.join(hou.expandString("$HIP"), "out", node.name() + '.json')
    # print(json.dumps(data, indent=4))
    print("    Save {} to {}".format(node.name(), filePath))
    with open(filePath, 'w') as f:
        f.write(json.dumps(data, ensure_ascii=False))

    print("    It takes {} s".format(time.time() - beg_time))

# def exportObjGeometry(node): end

## entry point
objs = hou.node("/obj")

for node in objs.recursiveGlob("*", hou.nodeTypeFilter.ObjGeometry):
    name = node.name()
    print("ObjGeometry: {}".format(name))

    if node.isObjectDisplayed() == False:
        print("    This isn't displayed.")
        continue

    exportObjGeometry(node)

print("--done--")
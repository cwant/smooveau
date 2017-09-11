import bpy, bmesh
from math import sqrt, sin, cos, atan2, pi
from mathutils import Matrix, Vector

class VertData:
  def __init__(self, v, bm, **kwargs):
    self.bvert = v
    self.input_vert = Vector(v.co)
    self.new_vert = None
    self.vert_index = v.index
    self.bm = bm
    self.input_edge_verts = []
    self.normal = None

  def add_edge_vert(self, edge, v):
    vec = Vector(v.co)
    edge_vert = { 'v': vec,
                  'i': v.index,
                  'e': edge,
                  'l': (self.input_vert - vec).magnitude }

    self.input_edge_verts.append(edge_vert)

  def vec_average(self, vecs):
    vave = Vector((0.0, 0.0, 0.0))
    for vec in vecs:
      vave += vec
    return vave * (1.0/len(vecs))

  def least_squares_normal(self, vave, diffs):
    # Method taken from here:
    # http://www.ilikebigbits.com/blog/2015/3/2/plane-from-points
    # Get normal vector for plane that best fits the vectors in diffs
    n = len(diffs)
    centroid = vave

    # Calc full 3x3 covariance matrix, excluding symmetries:
    (xx, xy, xz, yy, yz, zz) = [0.0]*6

    for p in diffs:
      r = p - centroid;
      xx += r[0] * r[0]
      xy += r[0] * r[1]
      xz += r[0] * r[2]
      yy += r[1] * r[1]
      yz += r[1] * r[2]
      zz += r[2] * r[2]

    det_x = yy*zz - yz*yz
    det_y = xx*zz - xz*xz
    det_z = xx*yy - xy*xy

    det_max = max([abs(det_x), abs(det_y), abs(det_z)])
    if (det_max < 0.00001):
      # Give up
      return vave

    if det_max == abs(det_x):
      a = (xz*yz - xy*zz) / det_x
      b = (xy*yz - xz*yy) / det_x
      normal = Vector((1.0, a, b))
    elif det_max == abs(det_y):
      a = (yz*xz - xy*zz) / det_y
      b = (xy*xz - yz*xx) / det_y
      normal = Vector((a, 1.0, b))
    else:
      a = (yz*xy - xz*yy) / det_z
      b = (xz*xy - yz*xx) / det_z
      normal = Vector((a, b, 1.0))
    if vave * normal < 0.0:
      return -normal
    return normal

  def calculate_shift(self):
    if len(self.input_edge_verts) == 0:
      return

    initial_length = sum([ev['l'] for ev in self.input_edge_verts])
    # Unit vectors pointing in direction of connected edges
    diffs = [(ev['v'] - self.input_vert) for ev in self.input_edge_verts]
    vave = self.vec_average(diffs)
    # Average point of verts at other end of edges
    ave = self.input_vert + vave
    ave_length = sum([(ev['v'] - ave).magnitude
                      for ev in self.input_edge_verts])

    # Things break if vave is really small (ambiguous curvature)
    if vave.magnitude > 0.0000001:
      if (len(self.input_edge_verts) < 3):
        # This is first guess at unit vector in pole direction
        self.normal = -vave.normalized()
      else:
        # This is second (better) guess
        self.normal = -self.least_squares_normal(vave, diffs).normalized()

    if ave_length < initial_length:
      ratio = sqrt((initial_length*initial_length) -
                   (ave_length * ave_length)) / len(self.input_edge_verts)
      ave += self.normal * ratio

    # Store for later (can't shift vert until other verts are calculated)
    self.new_vert = ave

  def shift_vertex(self):
    for i in range(3):
      self.bvert.co[i] = self.new_vert[i]

class Smooveau:
  def __init__(self, mesh, **kwargs):
    self.mesh = mesh

    self.kwargs = kwargs

    self.vert_data = []
    self.profile_connectors = []

  def modify_mesh(self):
    bm = bmesh.new()
    bm.from_mesh(self.mesh)

    self.create_vert_data(bm)
    self.add_edges_to_vert_data(bm)
    self.calculate_shift()
    self.shift_vertices()

    bm.to_mesh(self.mesh)
    bm.free()

  ### Below should be private

  def create_vert_data(self, bm):
    for vert in bm.verts:
      self.vert_data.append(VertData(vert, bm, **self.kwargs))

  def add_edges_to_vert_data(self, bm):
    for edge in bm.edges:
      for i in range(2):
        this_vert_index = edge.verts[i].index
        other_vert = edge.verts[(i + 1) % 2]
        self.vert_data[this_vert_index].add_edge_vert(edge, other_vert)

  def calculate_shift(self):
    for vert_data in self.vert_data:
      vert_data.calculate_shift()

  def shift_vertices(self):
    for vert_data in self.vert_data:
      vert_data.shift_vertex()

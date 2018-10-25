import bpy
import re
import bmesh
from mathutils import Vector, Matrix

# Copyright (C) 2018 Grzegorz Sobczak
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# grzeg.sobczak(at)gmail.com
# 
# MOL file structure based on http://www.daylight.com/meetings/mug05/Kappler/ctfile.pdf


def read_some_data(context, filepath, atom_dia, bond_dia):
    #print("Opening mol file...")
    f = open(filepath, 'r', encoding='utf-8')
    data = []
    for line in f:
        if line.strip():
            data.append(line)
    f.close()

    # would normally load the data here
    
    atomNum=int(data[1][0:3])
    bondNum=int(data[1][3:6])
    atomlist=int(data[1][6:9])
    chiral=[False, True][int(data[1][12:15])]
    version=data[1][33:36]
    
    scene = bpy.context.scene
    
    
    #ATOMS IMPORT
    
    atoms = []
    

    
    for i in range(2,2+atomNum):
        atom = {
            "x": 0.0,
            "y": 0.0,
            "z": 0.0,
            "name": "XXX",
            "charge": 0,
            "mesh": None,
            "object": None
        }
        atom["x"] = float(data[i][0:10])
        atom["y"] = float(data[i][10:20])
        atom["z"] = float(data[i][20:30])
        atom["name"]=data[i][31:34].strip()+str(i)
        atom["charge"]=int(data[i][36:39])
        atom["mesh"]=bpy.data.meshes.new(atom["name"]+"_mesh")
        atom["object"]=bpy.data.objects.new(atom["name"], atom["mesh"])
        scene.objects.link(atom["object"])
        scene.objects.active = atom["object"]
        atom["object"].select = True
        bm = bmesh.new()
        bmesh.ops.create_uvsphere(bm, u_segments=16, v_segments=16, diameter=atom_dia/2)
        bmesh.ops.translate(bm, vec=(atom["x"],atom["y"],atom["z"]), space=bpy.context.object.matrix_world, verts=bm.verts)
        bm.to_mesh(atom["mesh"])
        bm.free()
        bpy.ops.object.modifier_add(type='SUBSURF')
        bpy.ops.object.shade_smooth()
        #print(atom)
        
        atom["object"]["Element"] = re.findall("\D+", atom["name"])[0]
        atom["object"]["Charge"] = atom["charge"]
        
        atoms.append(atom)
    
    #BONDS IMPORT
    
    #print(atoms)
    
    bonds = []
    

    for i in range(2+atomNum, 2+atomNum+bondNum):
        bond = {
           "atom1": 0,
            "atom2": 0,
            "type": 1,
            "stereo":0,
            "topology": 0,
            "rcent":0,
            "mesh": None,
            "object": None
        }
        bond["atom1"]=int(data[i][0:3])
        bond["atom2"]=int(data[i][3:6])
        #print(atoms[bond["atom1"]-1]["name"])
        #print(atoms[bond["atom2"]-1]["name"])
        bond["type"]=int("0"+data[i][6:9].strip())
        bond["stereo"]=int("0"+data[i][9:12].strip())
        bond["topology"]=int("0"+data[i][15:18].strip())
        bond["rcent"]=int("0"+data[i][18:21].strip())
        bond["mesh"]=bpy.data.meshes.new("bond_"+atoms[bond["atom1"]-1]["name"]+"-"+atoms[bond["atom2"]-1]["name"]+"_mesh")
        bond["object"]=bpy.data.objects.new("bond_"+atoms[bond["atom1"]-1]["name"]+"-"+atoms[bond["atom2"]-1]["name"], bond["mesh"])
        scene.objects.link(bond["object"])
        scene.objects.active = bond["object"]
        bond["object"].select = True
        bm = bmesh.new()
        atom1_pos = Vector((atoms[bond["atom1"]-1]["x"], atoms[bond["atom1"]-1]["y"], atoms[bond["atom1"]-1]["z"]))
        atom2_pos = Vector((atoms[bond["atom2"]-1]["x"], atoms[bond["atom2"]-1]["y"], atoms[bond["atom2"]-1]["z"]))
        diff = atom2_pos-atom1_pos
        mean = (atom2_pos+atom1_pos)/2
        
        bmesh.ops.create_cone(bm, cap_ends=False, cap_tris=False, segments=16, diameter1=bond_dia/2, diameter2=bond_dia/2, depth=diff.length)
        bmesh.ops.translate(bm, vec=mean, space=bpy.context.object.matrix_world, verts=bm.verts)
        #cos=atom1_pos.dot(atom2_pos)/(atom1_pos.length*atom2_pos.length)
        #sin=atom1_pos.cross(atom2_pos)/(atom1_pos.length*atom2_pos.length)
        norm=Vector((0,0,1))
        rMatrix = Matrix.Rotation(norm.angle(diff), 4, norm.cross(diff))
        bmesh.ops.rotate(bm, cent=mean, matrix=rMatrix, verts=bm.verts)
        bm.to_mesh(bond["mesh"])
        bm.free
        
        bond["object"]["Atom 1"] = atoms[bond["atom1"]-1]["object"]
        bond["object"]["Atom 2"] = atoms[bond["atom2"]-1]["object"]
        bond["object"]["Order"] = bond["type"]
        
        bonds.append(bond)
        bond.clear()
    
    return {'FINISHED'}


# ImportHelper is a helper class, defines filename and
# invoke() function which calls the file selector.
from bpy_extras.io_utils import ImportHelper
from bpy.props import StringProperty, FloatProperty
from bpy.types import Operator


class ImportMolFile(Operator, ImportHelper):
    """This appears in the tooltip of the operator and in the generated docs"""
    bl_idname = "import_blenderchem_mol.data"
    bl_label = "Import MOL file"

    # ImportHelper mixin class uses this
    filename_ext = ".mol"

    filter_glob = StringProperty(
            default="*.mol",
            options={'HIDDEN'},
            maxlen=255,  # Max internal buffer length, longer would be clamped.
            )

    # List of operator properties, the attributes will be assigned
    # to the class instance from the operator settings before calling.
    atom_dia = FloatProperty(
            name="Atom diameter",
            description="Diameter of spheres representing atoms",
            min=0,
            soft_min=0,
            soft_max=2,
            precision=4,
            default=1.0,
            )

    bond_dia = FloatProperty(
            name="Bond thickness",
            description="Diameter of cylinders representing atoms",
            min=0,
            soft_min=0,
            soft_max=2,
            precision=4,
            default=0.4,
            )

    def execute(self, context):
        return read_some_data(context, self.filepath, self.atom_dia, self.bond_dia)


# Only needed if you want to add into a dynamic menu
def menu_func_import(self, context):
    self.layout.operator(ImportMolFile.bl_idname, text="Import MOL file")


def register():
    bpy.utils.register_class(ImportMolFile)
    bpy.types.INFO_MT_file_import.append(menu_func_import)


def unregister():
    bpy.utils.unregister_class(ImportMolFile)
    bpy.types.INFO_MT_file_import.remove(menu_func_import)


if __name__ == "__main__":
    register()

    # test call
    bpy.ops.import_blenderchem_mol.data('INVOKE_DEFAULT')

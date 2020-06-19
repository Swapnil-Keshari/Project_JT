# Project_JT

Added comments to the dataprocessing file
Added comments to the dataprocessing_1 file



Rewiring Network Algorithm

If mono exist
  -All mono to source
  If di exist
Check mono in di
If yes check tri exist
If tri exist check that di is in tri 
If yes the edge from tri to finish
If no: edge from di to tri
If no: then edge from di to tri
If mono is not in di then edge from mono to finish
If no di but tri exist
Check mono in tri
If exist edge from tri to finish
If not in tri then mono to finish
If no di tri exist
mono to finish
If no mono but di exist
All di to source
Check tri exist
Check di in tri
If yes then tri to finish
If no the di to finish
If tri doesn't exist then di to finish
If no mono and di exist
Add tri to source
Add edge tri to finish

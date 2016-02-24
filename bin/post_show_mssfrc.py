###import the simple module from the paraview
from paraview.simple import *


#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get active source.
h5PartReader32 = GetActiveSource()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1107, 789]

# get display properties
h5PartReader32Display = GetDisplayProperties(h5PartReader32, view=renderView1)

# hide color bar/color legend
h5PartReader32Display.SetScalarBarVisibility(renderView1, False)

for i in range(0, 32):
   name='H5PartReader'+str(i+1)
      
   # find source
   Render= FindSource(name)

   # set active source
   SetActiveSource(Render)

   # get active view
   renderView1 = GetActiveViewOrCreate('RenderView')

   # get display properties
   Disp = GetDisplayProperties(Render, view=renderView1)

   # set scalar coloring
   ColorBy(Disp, ('POINTS', 'mssfrc'))

   # rescale color and/or opacity maps used to include current data range
   Disp.RescaleTransferFunctionToDataRange(True)

   # show color bar/color legend
   Disp.SetScalarBarVisibility(renderView1, True)

   # get color transfer function/color map for 'Vx'
   mssfrcLUT = GetColorTransferFunction('mssfrc')

   # get opacity transfer function/opacity map for 'mssfrc'
   mssfrcPWF = GetOpacityTransferFunction('mssfrc')

   #

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get active source.
h5PartReader1 = GetActiveSource()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1107, 789]

# get display properties
h5PartReader1Display = GetDisplayProperties(h5PartReader1, view=renderView1)

# show color bar/color legend
h5PartReader1Display.SetScalarBarVisibility(renderView1, True)


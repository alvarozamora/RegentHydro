import "regent"

-- Helper modules to handle PNG files and command line arguments
local EdgeConfig = require("edge_config")
local coloring   = require("coloring_util")

-- Some C APIs
local c     = regentlib.c
local sqrt  = regentlib.sqrt(double)
local cmath = terralib.includec("math.h")
local PI = cmath.M_PI

-- Field space for Grid
fspace Grid2D
{
  --Grids
  U0  : double;   
  U1  : double;
  U2  : double;
  U3  : double; 

  --Fluxes
  F0  : double;
  F1  : double;
  F2  : double;
  F3  : double;
  G0  : double;
  G1  : double;
  G2  : double;
  G3  : double;

  --Updates
  X0  : double;
  X1  : double;
  X2  : double;
  X3  : double;

  Y0  : double;
  Y1  : double;
  Y2  : double;
  Y3  : double;
  

  c : double;
}

fspace fmax
{
  Max : double
}

fspace time
{
  Time : double;
  Cent : int32;
}

terra read_line(f : &c.FILE, my_struct_value : &double, size : int32, num : int32)
  return c.fread(my_struct_value, size, num, f)
end
terra dumphead(f : &c.FILE, val : &double)
  c.fwrite(val,8,16,f)
end

terra dump(f : &c.FILE, val : double)
  var a : double[1]
  a[0] = val
  c.fwrite(&a, 8, 1, f)
end


task factorize2d(parallelism : int) : int2d
  var limit = [int](cmath.sqrt([double](parallelism)))
  var size_x = 1
  var size_y = parallelism
  for i = 1, limit + 1 do
    if parallelism % i == 0 then
      size_x, size_y = i, parallelism / i
      if size_x > size_y then
        size_x, size_y = size_y, size_x
      end
    end
  end
  return int2d { size_x, size_y }
end

--
-- The 'initialize' task reads the image data from the file and initializes
-- the fields for later tasks. The upper left and lower right corners of the image
-- correspond to point {0, 0} and {width - 1, height - 1}, respectively.
--
task initialize(r_grid : region(ispace(int2d), Grid2D), size : int2d) -- TODO add file reader
where
  reads writes(r_grid)
do
  var KH : bool = true
  var load : bool = false

  --var delta : double = 0.02 -- produces nans
  --var delta : double = 0.000001 -- should recreate sharp interface 
  --var delta : double = 0.05 -- Tom's Paper
  var delta : double = 0.03 -- Testing Boundary
  if KH then  
    for e in r_grid do
      -- this is where I initialize values
      -- Kelvin-Helmholtz
      r_grid[e].U0 = 1 + 1/(1+cmath.exp(-2*(e.y*1.0/size.y-0.25)/delta))/(1+cmath.exp(-2*(0.75-e.y*1.0/size.y)/delta))
      r_grid[e].U1 = (-0.5 + 1/(1+cmath.exp(-2*(e.y*1.0/size.y-0.25)/delta))/(1+cmath.exp(-2*(0.75-e.y*1.0/size.y)/delta)))*r_grid[e].U0
      r_grid[e].U2 = 0.05*cmath.sin(2*PI*e.x/size.x)*r_grid[e].U0
      r_grid[e].U3 = 2.5/(5./3-1) + 1./2*(r_grid[e].U1*r_grid[e].U1 + r_grid[e].U2*r_grid[e].U2)/r_grid[e].U0 --Pressure Eq
      
      --c.printf("U = {%f, %f, %f, %f}\n", r_grid[e].U0, r_grid[e].U1, r_grid[e].U2, r_grid[e].U3)
      
      --if e.y >= 3*size.y/4 or e.y <= size.y/4 then
      --  r_grid[e].U0 = 1.
      --  r_grid[e].U1 = -0.5*r_grid[e].U0
      --  r_grid[e].U2 = 0.04*cmath.sin(2*PI*e.x/size.x)*r_grid[e].U0
      --  --r_grid[e].U3 = (r_grid[e].U1*r_grid[e].U1 + r_grid[e].U2*r_grid[e].U2)/r_grid[e].U0
      --  r_grid[e].U3 = 2.5/(5./3-1) + 1./2*(r_grid[e].U1*r_grid[e].U1 + r_grid[e].U2*r_grid[e].U2)/r_grid[e].U0 --Pressure Eq
      --else 
      --  r_grid[e].U0 = 2.
      --  r_grid[e].U1 = 0.5*r_grid[e].U0
      --  r_grid[e].U2 = 0.04*cmath.sin(2*PI*e.x/size.x)*r_grid[e].U0 
      --  --r_grid[e].U3 = (r_grid[e].U1*r_grid[e].U1 + r_grid[e].U2*r_grid[e].U2)/r_grid[e].U0
      --  r_grid[e].U3 = 2.5/(5./3-1) + 1./2*(r_grid[e].U1*r_grid[e].U1 + r_grid[e].U2*r_grid[e].U2)/r_grid[e].U0 --Pressure Eq
      --end      



      --if e.y >= size.y/2 then
      --  r_grid[e].U0 = 1.0-0.225
      --  r_grid[e].U1 = -0.225*r_grid[e].U0
      --  r_grid[e].U2 = 0.025*cmath.sin(2*PI*e.x/size.x)*r_grid[e].U0
      --  --r_grid[e].U3 = (r_grid[e].U1*r_grid[e].U1 + r_grid[e].U2*r_grid[e].U2)/r_grid[e].U0
      --  r_grid[e].U3 = 0.4/(5./3-1) + 1./2*(r_grid[e].U1*r_grid[e].U1 + r_grid[e].U2*r_grid[e].U2)/r_grid[e].U0 --Pressure Eq
      --else 
      --  r_grid[e].U0 = 1.0+0.225
      --  r_grid[e].U1 = 0.225*r_grid[e].U0
      --  r_grid[e].U2 = 0.025*cmath.sin(2*PI*e.x/size.x)*r_grid[e].U0 
      --  --r_grid[e].U3 = (r_grid[e].U1*r_grid[e].U1 + r_grid[e].U2*r_grid[e].U2)/r_grid[e].U0
      --  r_grid[e].U3 = 0.4/(5./3-1) + 1./2*(r_grid[e].U1*r_grid[e].U1 + r_grid[e].U2*r_grid[e].U2)/r_grid[e].U0 --Pressure Eq
      --end
    end
  end
  



  if load then
    var f = c.fopen('data0.npy', 'rb')
    var head : double[16]
    read_line(f, head, 8, 16)
    for e in r_grid do
      var data : double[1]
      read_line(f, data, 8, 1)
      r_grid[e].U0 = data[0]
    end
    for e in r_grid do
      var data : double[1] 
      read_line(f, data, 8, 1)
      r_grid[e].U1 = (data[0]-1)*r_grid[e].U0
    end
    for e in r_grid do
      var data : double[1]
      read_line(f, data, 8, 1)
      r_grid[e].U2 = (data[0]-1)*r_grid[e].U0
    end  
    for e in r_grid do
      var data : double[1]
      read_line(f, data, 8, 1)
      r_grid[e].U3 = data[0]
    end  
  end




  return 1
end

terra terraSign(x : double)
  if x >= 0 then
    return 1
  else
    return -1 
  end
end

terra terraAbs(x : double)
  if x >= 0 then
    return x
  else
    return -x
  end
end


terra terraAbsMin(x : double, y : double, z : double)
  if terraAbs(x) <= terraAbs(y) and terraAbs(x) <= terraAbs(z) then
    return terraAbs(x)
  elseif terraAbs(y) <= terraAbs(x) and terraAbs(y) <= terraAbs(z) then
    return terraAbs(y)
  elseif terraAbs(z) <= terraAbs(x) and terraAbs(z) <= terraAbs(y) then
    return terraAbs(z)
  end 
end

terra terraminmod(x : double, y: double, z: double)
  return 1/4.*terraAbs(terraSign(x) + terraSign(y))*(terraSign(x) + terraSign(z))*terraAbsMin(x,y,z)
end

terra terraPLM(i1 : double, i2 : double, i3 : double, right : bool)
  var t : double = 1.0 --theta
  if right then
    return i2 - 0.5*terraminmod(t*(i2-i1), 0.5*(i3-i1), t*(i3-i2))
  else 
    return i2 + 0.5*terraminmod(t*(i2-i1), 0.5*(i3-i1), t*(i3-i2))
  end
end

task UpdateFlux(r_grid : region(ispace(int2d), Grid2D),
                g : double)
where
  reads(r_grid.{U0,U1,U2,U3,F1,F2,F3}),
  writes(r_grid.{F0,F1,F2,F3,G0,G1,G2,G3,c})
do
  var ts_start = c.legion_get_current_time_in_micros()
  for e in r_grid do    
    --Update Fluxes
    r_grid[e].F0 = r_grid[e].U1
    r_grid[e].F1 = r_grid[e].U1 * r_grid[e].U1 / r_grid[e].U0 + (g-1)*(r_grid[e].U3 - 1./2*(r_grid[e].U1*r_grid[e].U1 + r_grid[e].U2*r_grid[e].U2)/r_grid[e].U0)
    r_grid[e].F2 = r_grid[e].U1 * r_grid[e].U2 / r_grid[e].U0
    r_grid[e].F3 = (r_grid[e].U3 + (g-1)*(r_grid[e].U3 - 1./2*(r_grid[e].U1*r_grid[e].U1 + r_grid[e].U2*r_grid[e].U2)/r_grid[e].U0))*r_grid[e].U1/r_grid[e].U0 
    r_grid[e].G0 = r_grid[e].U2
    r_grid[e].G1 = r_grid[e].F2
    r_grid[e].G2 = r_grid[e].F1 - (r_grid[e].U1*r_grid[e].U1 - r_grid[e].U2*r_grid[e].U2)/r_grid[e].U0
    r_grid[e].G3 = r_grid[e].F3/r_grid[e].U1*r_grid[e].U2

    --Update Speed
    r_grid[e].c = cmath.sqrt(g*(g-1)*(r_grid[e].U3 - 1./2*(r_grid[e].U1*r_grid[e].U1 + r_grid[e].U2*r_grid[e].U2)/r_grid[e].U0)/r_grid[e].U0)

  end
  var ts_end = c.legion_get_current_time_in_micros()
  --c.printf("UpdateFlux took %.5f sec.\n", (ts_end - ts_start) * 1e-6)
end

task HLL_x(r_grid : region(ispace(int2d), Grid2D),
           r_max : region(ispace(int2d), fmax),
           g : double, dx : double, size : int2d)
where
  reads atomic(r_grid.{U0,U1,U2,U3,F0,F1,F2,F3,c}),
  reads writes(r_grid.{X0,X1,X2,X3}),
  writes(r_max)
do
  var ts_start = c.legion_get_current_time_in_micros()
  var amax : double = 0
  for e in r_grid do
      --Calculate Eigenvalues for F
      var ap : double = 0
      var am : double = 0

      var left = e-{1,0}
      if left.x < 0 then
        left = {left.x+size.x, left.y}
      end
      var right = e+{1,0}
      if right.x >= size.x then
        right = {right.x-size.x, right.y}
      end
      var right2 = e+{2,0}
      if right2.x >= size.x then
        right2 = {right2.x-size.x, right2.y}
      end
      -- if c + vx > 0 set ap to c+vx
      if (r_grid[e].c + r_grid[e].U1/r_grid[e].U0) > 0 then
        ap = r_grid[e].c + r_grid[e].U1/r_grid[e].U0
      end

      --if c + vx  on the right cell is > left, set ap to c+vx
      if (r_grid[right].c + r_grid[right].U1/r_grid[right].U0) > ap then
        ap = r_grid[right].c + r_grid[right].U1/r_grid[right].U0 
      end
      if ap > amax then
        amax = ap
      end

      -- if c - vx > 0 set am to c-vx
      if (r_grid[e].c - r_grid[e].U1/r_grid[e].U0) > 0 then
        am = r_grid[e].c - r_grid[e].U1/r_grid[e].U0
      end
      -- if c - vx on the right cell is > left, set am to c-vx
      if (r_grid[right].c - r_grid[right].U1/r_grid[right].U0) > am then
        am = r_grid[right].c - r_grid[right].U1/r_grid[right].U0
      end
      if am > amax then
        amax = am
      end
      --Right HLL flux, Compute PLM values
      var FL0 : double = terraPLM(r_grid[left].F0,r_grid[e].F0,r_grid[right].F0,false)
      var FL1 : double = terraPLM(r_grid[left].F1,r_grid[e].F1,r_grid[right].F1,false)
      var FL2 : double = terraPLM(r_grid[left].F2,r_grid[e].F2,r_grid[right].F2,false)
      var FL3 : double = terraPLM(r_grid[left].F3,r_grid[e].F3,r_grid[right].F3,false)

      var FR0 : double = terraPLM(r_grid[e].F0,r_grid[right].F0,r_grid[right2].F0,true)
      var FR1 : double = terraPLM(r_grid[e].F1,r_grid[right].F1,r_grid[right2].F1,true)
      var FR2 : double = terraPLM(r_grid[e].F2,r_grid[right].F2,r_grid[right2].F2,true)
      var FR3 : double = terraPLM(r_grid[e].F3,r_grid[right].F3,r_grid[right2].F3,true)

      var UL0 : double = terraPLM(r_grid[left].U0,r_grid[e].U0,r_grid[right].U0,false)
      var UL1 : double = terraPLM(r_grid[left].U1,r_grid[e].U1,r_grid[right].U1,false)
      var UL2 : double = terraPLM(r_grid[left].U2,r_grid[e].U2,r_grid[right].U2,false)
      var UL3 : double = terraPLM(r_grid[left].U3,r_grid[e].U3,r_grid[right].U3,false)

      var UR0 : double = terraPLM(r_grid[e].U0,r_grid[right].U0,r_grid[right2].U0,true)
      var UR1 : double = terraPLM(r_grid[e].U1,r_grid[right].U1,r_grid[right2].U1,true)
      var UR2 : double = terraPLM(r_grid[e].U2,r_grid[right].U2,r_grid[right2].U2,true)
      var UR3 : double = terraPLM(r_grid[e].U3,r_grid[right].U3,r_grid[right2].U3,true)

      var HLL0 : double = (ap*FL0 + am*FR0- ap*am*(UR0 - UL0))/(ap + am)
      var HLL1 : double = (ap*FL1 + am*FR1- ap*am*(UR1 - UL1))/(ap + am)
      var HLL2 : double = (ap*FL2 + am*FR2- ap*am*(UR2 - UL2))/(ap + am)
      var HLL3 : double = (ap*FL3 + am*FR3- ap*am*(UR3 - UL3))/(ap + am)

 
      --Reduce L
      r_grid[e].X0 = r_grid[e].X0-HLL0/dx
      r_grid[right].X0 = r_grid[right].X0+HLL0/dx
      r_grid[e].X1 = r_grid[e].X1-HLL1/dx
      r_grid[right].X1 = r_grid[right].X1+HLL1/dx
      r_grid[e].X2 = r_grid[e].X2-HLL2/dx
      r_grid[right].X2 = r_grid[right].X2+HLL2/dx
      r_grid[e].X3 = r_grid[e].X3-HLL3/dx
      r_grid[right].X3 = r_grid[right].X3+HLL3/dx

  end
  var ts_end = c.legion_get_current_time_in_micros()
  --c.printf("HLL_x took %.5f sec.\n", (ts_end - ts_start) * 1e-6)
  __fence(__execution, __block)
  for e in r_max do
    r_max[e].Max = amax
  end
end


task HLL_y(r_grid : region(ispace(int2d), Grid2D),
           r_max : region(ispace(int2d), fmax),
           g : double, dy : double, size : int2d)
where
  reads atomic(r_grid.{U0,U1,U2,U3,G0,G1,G2,G3,c}),
  reads writes(r_grid.{Y0,Y1,Y2,Y3}),
  writes(r_max)
do
  var ts_start = c.legion_get_current_time_in_micros()
  var amax : double = 0
    for e in r_grid do
      --Calculate Eigenvalues for G
      var ap : double = 0
      var am : double = 0

      var left = e-{0,1}
      if left.y < 0 then
        left = {left.x, left.y+size.y}
      end
      var right = e+{0,1}
      if right.y >= size.y then
        right = {right.x, right.y-size.y}
      end
      var right2 = e+{0,2}
      if right2.y >= size.y then
        right2 = {right2.x, right2.y-size.y}
      end


      -- if c + vy > 0 set ap to c+vy
      if (r_grid[e].c + r_grid[e].U2/r_grid[e].U0) > 0 then
        ap = r_grid[e].c + r_grid[e].U2/r_grid[e].U0
      end

      --if c + vy > 0 on the right cell, set ap to c+vy
      if (r_grid[right].c + r_grid[right].U2/r_grid[right].U0) > ap then
        ap = r_grid[right].c + r_grid[right].U2/r_grid[right].U0
      end
      if ap > amax then
        amax = ap
      end
      -- if c - vy > 0 set am to c-vy
      if (r_grid[e].c - r_grid[e].U2/r_grid[e].U0) > 0 then
        am = r_grid[e].c - r_grid[e].U2/r_grid[e].U0
      end
      -- if c - vy > 0 on the right cell, set am to c+vy
      if (r_grid[right].c - r_grid[right].U2/r_grid[right].U0) > am then
        am = r_grid[right].c - r_grid[right].U2/r_grid[right].U0
      end
      if am > amax then
        amax = am
      end

      --Right HLL flux, First the PLM values
      var GL0 : double = terraPLM(r_grid[left].G0,r_grid[e].G0,r_grid[right].G0,false)
      var GL1 : double = terraPLM(r_grid[left].G1,r_grid[e].G1,r_grid[right].G1,false)
      var GL2 : double = terraPLM(r_grid[left].G2,r_grid[e].G2,r_grid[right].G2,false)
      var GL3 : double = terraPLM(r_grid[left].G3,r_grid[e].G3,r_grid[right].G3,false)

      var GR0 : double = terraPLM(r_grid[e].G0,r_grid[right].G0,r_grid[right2].G0,true)
      var GR1 : double = terraPLM(r_grid[e].G1,r_grid[right].G1,r_grid[right2].G1,true)
      var GR2 : double = terraPLM(r_grid[e].G2,r_grid[right].G2,r_grid[right2].G2,true)
      var GR3 : double = terraPLM(r_grid[e].G3,r_grid[right].G3,r_grid[right2].G3,true)

      var UL0 : double = terraPLM(r_grid[left].U0,r_grid[e].U0,r_grid[right].U0,false)
      var UL1 : double = terraPLM(r_grid[left].U1,r_grid[e].U1,r_grid[right].U1,false)
      var UL2 : double = terraPLM(r_grid[left].U2,r_grid[e].U2,r_grid[right].U2,false)
      var UL3 : double = terraPLM(r_grid[left].U3,r_grid[e].U3,r_grid[right].U3,false)

      var UR0 : double = terraPLM(r_grid[e].U0,r_grid[right].U0,r_grid[right2].U0,true)
      var UR1 : double = terraPLM(r_grid[e].U1,r_grid[right].U1,r_grid[right2].U1,true)
      var UR2 : double = terraPLM(r_grid[e].U2,r_grid[right].U2,r_grid[right2].U2,true)
      var UR3 : double = terraPLM(r_grid[e].U3,r_grid[right].U3,r_grid[right2].U3,true)

      var HLL0 : double = (ap*GL0 + am*GR0- ap*am*(UR0 - UL0))/(ap + am)
      var HLL1 : double = (ap*GL1 + am*GR1- ap*am*(UR1 - UL1))/(ap + am)
      var HLL2 : double = (ap*GL2 + am*GR2- ap*am*(UR2 - UL2))/(ap + am)
      var HLL3 : double = (ap*GL3 + am*GR3- ap*am*(UR3 - UL3))/(ap + am)

      --Reduce L
      r_grid[e].Y0 = r_grid[e].Y0-HLL0/dy
      r_grid[right].Y0 = r_grid[right].Y0+HLL0/dy
      r_grid[e].Y1 = r_grid[e].Y1-HLL1/dy
      r_grid[right].Y1 = r_grid[right].Y1+HLL1/dy
      r_grid[e].Y2 = r_grid[e].Y2-HLL2/dy
      r_grid[right].Y2 = r_grid[right].Y2+HLL2/dy
      r_grid[e].Y3 = r_grid[e].Y3-HLL3/dy
      r_grid[right].Y3 = r_grid[right].Y3+HLL3/dy
      

  end
  var ts_end = c.legion_get_current_time_in_micros()
  --c.printf("HLL_y took %.5f sec.\n", (ts_end - ts_start) * 1e-6)
  __fence(__execution, __block)  
  for e in r_max do
    r_max[e].Max = amax
  end
end

task MaxAlpha(r_xmax : region(ispace(int2d), fmax), 
              r_ymax : region(ispace(int2d), fmax))
where
  reads(r_xmax),
  reads(r_ymax)
do
  var MAX : double = 0
  for e in r_xmax do
    if r_xmax[e].Max > MAX then
      MAX = r_xmax[e].Max
    end
  end

  for e in r_ymax do  
    if r_ymax[e].Max > MAX then
      MAX = r_ymax[e].Max
    end
  end
  
  return MAX
end

task TimeStep(r_grid : region(ispace(int2d),Grid2D), d : double, MAX : double, r_time : region(ispace(int2d),time), endtime: double)
where
  reads(r_grid.{X0,X1,X2,X3,Y0,Y1,Y2,Y3}),
  reads writes(r_time),
  reads writes(r_grid.{U0,U1,U2,U3})
do
  var ts_start = c.legion_get_current_time_in_micros()

  --Calculates dt --time evolution
  var cf : double = 0.04 -- Courant Factor
  var dt : double = cf*d/MAX
  var frames : int32 = 300
  for e in r_time do
    if r_time[e].Time + dt > r_time[e].Cent*endtime/frames then
      dt = r_time[e].Cent*endtime/frames - r_time[e].Time
      r_time[e].Cent = r_time[e].Cent + 1
    end
    if r_time[e].Time + dt > endtime then
      dt = endtime - r_time[e].Time
    end

    r_time[e].Time += dt
  end
  
  for e in r_grid do
    r_grid[e].U0 = r_grid[e].U0 + dt*(r_grid[e].X0 + r_grid[e].Y0)
    r_grid[e].U1 = r_grid[e].U1 + dt*(r_grid[e].X1 + r_grid[e].Y1)
    r_grid[e].U2 = r_grid[e].U2 + dt*(r_grid[e].X2 + r_grid[e].Y2)
    r_grid[e].U3 = r_grid[e].U3 + dt*(r_grid[e].X3 + r_grid[e].Y3)
  end
  var ts_end = c.legion_get_current_time_in_micros()
  --c.printf("TimeStep took %.5f sec.\n", (ts_end - ts_start) * 1e-6)
end


terra wait_for(x : int) return 1 end
task block_task(r_image : region(ispace(int2d), Grid2D))
where
  reads writes(r_image)
do
  return 1
end


task Dump(r_grid : region(ispace(int2d), Grid2D), iter : int32)
where 
  reads writes(r_grid)
do
  var filename : int8[1000]
  c.sprintf([&int8](filename), 'Grid%d',iter)
  var g = c.fopen(filename,'wb')

  for e in r_grid do
    dump(g,r_grid[e].U0)
  end
  __fence(__execution, __block)
  for e in r_grid do
    dump(g,r_grid[e].U1)
  end
  __fence(__execution, __block)
  for e in r_grid do
    dump(g,r_grid[e].U2)
  end
  __fence(__execution, __block)
  for e in r_grid do
    dump(g,r_grid[e].U3)
  end
  
  c.fclose(g)
end


task toplevel()
  var config : EdgeConfig
  config:initialize_from_command()

  -- Create a logical region for Grid
  --var size_image = png.get_image_size(config.filename_image)
  var size : int2d = {2048,2048}
  var r_grid = region(ispace(int2d, size), Grid2D)
  var p_x = partition(equal, r_grid, ispace(int2d, {1, config.parallelism}))
  var p_y = partition(equal, r_grid, ispace(int2d, {config.parallelism, 1}))
  
  -- Create a halo partition for ghost access, will need for GP-based stencil
  --var c_halo = coloring.create()
  --for color in p_private_colors do
  --  var bounds = p_private[color].bounds
  --  var halo_bounds : rect2d = {bounds.lo - {1,1}, bounds.hi + {2,2}}-- Bounds are set according to spatial reconstruction scheme -- PLM is -1, +2
  --  coloring.color_domain(c_halo, color, halo_bounds)
  --end 
  
  --Create an aliased partition of region 'r_fake' using coloring 'c_halo':
  --var p_halo = partition(aliased, r_fake, c_halo, p_private_colors)
  --coloring.destroy(c_halo)
  
  --Create region to gather maximum eigenvalues
  var r_xmax = region(ispace(int2d, {1, config.parallelism}), fmax)  
  var r_ymax = region(ispace(int2d, {config.parallelism, 1}), fmax)  
  var p_xmax = partition(equal, r_xmax, p_x.colors)
  var p_ymax = partition(equal, r_ymax, p_y.colors)
  
  --Create a region to store dt
  var r_time = region(ispace(int2d, {config.parallelism, 1}), time)
  --var r_time = region(ispace(int2d, {1, config.parallelism}), time)
  var p_time = partition(equal, r_time, p_y.colors)

  var g : double = 5./3 -- Adiabatic Index
  var lx : uint32 = size.x
  var ly : uint32 = size.y
  var dx : double = 1./lx
  var dy : double = 1./ly
  var d : double = dx
  if dy < dx then
     d = dy
  end

  --var token = initialize(r_grid)
  --Use this scheme for initializing with KHI etc
  var token : int32 = 1
  for color in p_x.colors do
    initialize(p_x[color], size)
    token += block_task(p_x[color])
  end

  Dump(r_grid, 0)
  
  wait_for(token)
  var TS_start = c.legion_get_current_time_in_micros()
  
  var IO : bool = true 
  var SimTime : double = 0
  var endtime : double = 2.0
  var iter : int32 = 0
  var cent : int32 = 1
  var maxiter : int32 = 10000000
  var MAX : double = 0
  fill(r_time.Time, SimTime)  
  fill(r_time.Cent, cent)  
  while SimTime < endtime do -- and iter < maxiter do
  iter += 1

  --Update Flux, Sound Speed, and Zero L's
  for color in p_x.colors do
    UpdateFlux(p_x[color], g)
  end

  --Compute HLL flux and reduce L
  for color in p_x.colors do 
    fill((p_x[color]).X0, 0)
    fill((p_x[color]).X1, 0)
    fill((p_x[color]).X2, 0)
    fill((p_x[color]).X3, 0)
  end 
  for color in p_y.colors do 
    fill((p_y[color]).Y0, 0)
    fill((p_y[color]).Y1, 0)
    fill((p_y[color]).Y2, 0)
    fill((p_y[color]).Y3, 0)
  end 

  for color in p_x.colors do 
     HLL_x(p_x[color], p_xmax[color], g, dx,size)
  end 
  for color in p_y.colors do
     HLL_y(p_y[color], p_ymax[color], g, dy,size)
  end
  
  if (iter - 1)%500 == 0 then
    MAX = MaxAlpha(r_xmax,r_ymax)
  end

  var OldCent : int32 = r_time[{0,0}].Cent
  for color in p_y.colors do --need to adjust inputs
    TimeStep(p_y[color], d, MAX, p_time[color], endtime)
  end
  var LastTime : double = SimTime 
  SimTime = r_time[{0,0}].Time

  if r_time[{0,0}].Cent > OldCent then
    Dump(r_grid, OldCent)
  end 
  c.printf("Sim Time = %f, Last step = %f \n", SimTime, SimTime - LastTime)  

  end --this one ends the 'while T < endtime' loop
  
  for color in p_y.colors do
    token += block_task(p_y[color])
  end
  wait_for(token)
  var TS_end = c.legion_get_current_time_in_micros()
  c.printf("Total time: %.6f sec. Total Iterations: %d\n", (TS_end - TS_start) * 1e-6, iter)
  c.printf("Sim Time = %f\n", SimTime)
  --TODO Save Grid
  --var filename = "density"
  --attach(hdf5, r_fake.U0, filename, regentlib.file_read_write)  
  --release(r_fake.U0)
  --detach(hdf5, r_fake.U0)
  --print_density(r_grid[0])


  Dump(r_grid, iter)
end

regentlib.start(toplevel)

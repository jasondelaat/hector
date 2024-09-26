do
   -- Simple Types -------------------------------------------------------------
   local function create(meta)
      if not meta.__index then
         meta.__index = meta
      end
      return setmetatable({}, meta)
   end

   local function newtype()
      return setmetatable(
         {},
         {
            __call=function(self, ...)
               local obj = create(self)
               if obj.init then
                  obj:init(...)
               end
               return obj
            end
         }
      )
   end
   -- END Simple Types ---------------------------------------------------------










   -- Utilities --------------------------------------------------------
   local function contains(lst, val)
      for _,v in ipairs(lst) do
         if v == val then
            return true
         end
      end
   end

   local function get(tbl, index, default)
      local value = tbl[index]
      if not value then
         return default
      end
      return value
   end

   local function map(f, lst)
      local r = {}
      for i,v in ipairs(lst) do
         r[i] = f(v)
      end
      return r
   end

   local function map_keys(f, tbl)
      local r = {}
      for k,v in pairs(tbl) do
         r[k] = f(v, k)
      end
      return r
   end

   local function round(n)
      local i = math.floor(n)
      if n-i >= 0.5 then
         return i + 1
      else
         return i
      end
   end

   -- Square roots. These get used a bunch so just calculate them once.
   local sqrt2 = math.sqrt(2)
   local sqrt3 = math.sqrt(3)

   -- END Utilities ----------------------------------------------------










   -- Multivectors ---------------------------------------------------------

   local multivector -- forward declaration
   
   function rotor(a, b, c, d)
      return multivector{scalar=a, xy=b, xz=c, yz=d}
   end

   function vector(x, y, z)
      return multivector{x=x, y=y, z=z}
   end

   local components = {
      'scalar',
      'x', 'y', 'z',
      'xy', 'xz', 'yz',
      'xyz'
   }

   local mult_table = {xyz={xyz={-1, "scalar"},xz={1, "y"},scalar={1, "xyz"},yz={-1, "x"},y={-1, "xz"},x={1, "yz"},z={1, "xy"},xy={-1, "z"}},xz={xyz={1, "y"},xz={-1, "scalar"},scalar={1, "xz"},yz={-1, "xy"},y={-1, "xyz"},x={-1, "z"},z={1, "x"},xy={1, "yz"}},scalar={xyz={1, "xyz"},xz={1, "xz"},scalar={1, "scalar"},yz={1, "yz"},y={1, "y"},x={1, "x"},z={1, "z"},xy={1, "xy"}},yz={xyz={-1, "x"},xz={1, "xy"},scalar={1, "yz"},yz={-1, "scalar"},y={-1, "z"},x={1, "xyz"},z={1, "y"},xy={-1, "xz"}},y={xyz={-1, "xz"},xz={-1, "xyz"},scalar={1, "y"},yz={1, "z"},y={1, "scalar"},x={-1, "xy"},z={1, "yz"},xy={-1, "x"}},x={xyz={1, "yz"},xz={1, "z"},scalar={1, "x"},yz={1, "xyz"},y={1, "xy"},x={1, "scalar"},z={1, "xz"},xy={1, "y"}},z={xyz={1, "xy"},xz={-1, "x"},scalar={1, "z"},yz={-1, "y"},y={-1, "yz"},x={-1, "xz"},z={1, "scalar"},xy={1, "xyz"}},xy={xyz={-1, "z"},xz={-1, "yz"},scalar={1, "xy"},yz={1, "xz"},y={1, "x"},x={-1, "y"},z={1, "xyz"},xy={-1, "scalar"}}}

   local grades = {scalar=0,x=1,y=1,z=1,xy=2,xz=2,yz=2,xyz=3}

   multivector = newtype()
   function multivector:init(cs)
      cs = cs or {}
      for _,c in ipairs(components) do
         local val = cs[c] or 0
         self[c] = math.abs(val) > 0.000000001 and val or nil
      end
   end

   function multivector:abs()
      return multivector(map_keys(math.abs, self))
   end

   function multivector:floor()
      return multivector(map_keys(math.floor, self))
   end

   function multivector.hex_cross(a, b)
      local p = hex_to_standard(a)
      local q = hex_to_standard(b)
      return sqrt3*(p^q) / 2
   end

   function multivector.hex_dot(a, b)
      local p = hex_to_standard(a)
      local q = hex_to_standard(b)
      return (p..q) / 2
   end

   function multivector:hex_len()
      return (
         self:taxi_len() + math.abs(hector.constraint_factor*self.x - self.y)
      )/2
   end

   function multivector:magnitude()
      return math.sqrt(self:magnitude2())
   end

   function multivector:magnitude2()
      return (self*self:reverse()).scalar
   end

   function multivector:neighbours()
      local std = hex_to_standard(self)
      local ns = {}
      for _,v in ipairs({P_AXIS, Q_AXIS, R_AXIS}) do
         table.insert(ns, std + v)
         table.insert(ns, std - v)
      end
      return map(standard_to_hex, ns)
   end
   multivector.neighbors = multivector.neighbours

   function multivector:reflect(v)
      v = v / sqrt2
      local std = hex_to_standard(self)
      return standard_to_hex(v*std*v)
   end

   function multivector:reverse()
      return multivector(map_keys(function(v, k) return grades[k] > 1 and -v or v end, self))
   end

   function multivector:rotate(n, center)
      n = n or 1
      center = center or vector(0, 0)
      local v = hex_to_standard(self - center)
      local r = hector.rotors.sixty
      if n < 0 then
         r = r:reverse()
      end
      for i=1,math.abs(n)-1 do
         r = r * r
      end
      return standard_to_hex(r:reverse()*v*r) + center
   end

   function multivector:round()
      local x = round(self.x)
      local y = round(self.y)
      local z = round(self.z)
      local dx = math.abs(self.x - x)
      local dy = math.abs(self.y - y)
      local dz = math.abs(self.z - z)
      if dx > dy and dx > dz then
         return vector(-y-z, y, z)
      elseif dy > dx and dy > dz then
         return vector(x, -x-z, z)
      else
         return vector(x, y, -x-y)
      end
   end

   function multivector:screen_len()
      local p = self:to_screen()
      return math.sqrt(p.x*p.x + p.y*p.y)
   end

   function multivector:show()
      local s = ''
      for _,c in ipairs(components) do
         local v = self[c]
         if v ~= 0 then
            local sign = v > 0 and '+' or '-'
            local comp = c == 'scalar' and '' or c
            s = s..string.format(' %s %f%s', sign, math.abs(v), comp)
         end
      end
      return s
   end

   function multivector:taxi_len()
      local temp = self:abs()
      return temp.x + temp.y + temp.z
   end

   function multivector:to_hex()
      local r = hector.rotors[hector.style]
      local frac = vector(self.x / hector.scale_x, self.y / hector.scale_y, 0)
      local std = r*frac*r:reverse()
      return standard_to_hex(std:round())
   end

   function multivector:to_screen()
      local std = hex_to_standard(self)
      local r = hector.rotors[hector.style]
      local rot = r:reverse()*std*r
      return vector(
         hector.scale_x*rot.x,
         hector.scale_y*rot.y
      ):floor()
   end

   function multivector.__add(a, b)
      local cs = {}
      for _,c in ipairs(components) do
         cs[c] = a[c] + b[c]
      end
      return multivector(cs)
   end

   -- Not a full GA dot product just a normal 3D vector dot product.
   function multivector.__concat(a, b)
      return a.x*b.x + a.y*b.y + a.z*b.z
   end

   function multivector:__div(n)
      local cs = {}
      for _,c in ipairs(components) do
         cs[c] = self[c] / n
      end
      return multivector(cs)
   end

   function multivector.__eq(a, b)
      for _,c in ipairs(components) do
         if a[c] ~= b[c] then
            return false
         end
      end
      return true
   end

   function multivector:__index(key)
      local val = rawget(self, key)
      if contains(components, key) and not val then
         return 0
      else
         return val or multivector[key]
      end
   end

   function multivector.__mul(a, b)
      if type(a) == 'number' then
         a = multivector{scalar=a}
      end
      local results = {}
      for c1,v1 in pairs(a) do
         local cs = {}
         for c2,v2 in pairs(b) do
            local sign, comp = table.unpack(mult_table[c1][c2])
            cs[comp] = sign * v1 * v2
         end
         table.insert(results, multivector(cs))
      end
      if #results == 0 then
         return multivector()
      else
         local r = results[1]
         for i=2,#results do
            r = r + results[i]
         end
         return r
      end
   end

   function multivector.__pow(a, b)
      return a.y*b.x - a.x*b.y
   end

   function multivector.__sub(a, b)
      local cs = {}
      for _,c in ipairs(components) do
         cs[c] = a[c] - b[c]
      end
      return multivector(cs)
   end

   function multivector:__tostring()
      return string.format('(%.0f, %.0f)', self.x, self.y)
   end

   function multivector:__unm()
      local cs = {}
      for _,c in ipairs(components) do
         cs[c] = -self[c]
      end
      return multivector(cs)
   end


   -- END Multivectors -----------------------------------------------------










   -- Hector Initialization ----------------------------------------------------

   -- Global Constants
   P_AXIS = vector(1, 0, -1)
   Q_AXIS = vector(-1, 1, 0)
   R_AXIS = vector(0, -1, 1)

   FLAT = 'flat'
   POINTY = 'pointy'

   -- Converting to/from standard coordinates. (AKA cube coordinates)
   function hex_to_standard(hex)
      return hex.x * hector.x_axis + hex.y * hector.y_axis
   end

   function standard_to_hex(std)
      return vector(
         std..hector.x_inv,
         std..hector.y_inv
      )
   end

   hector = newtype()
   function hector.init(options)
      options = options or {}
      hector.style = get(options, 'style', FLAT)
      hector.x_axis = get(options, 'x_axis', P_AXIS)
      hector.y_axis = get(options, 'y_axis', -R_AXIS)
      hector.major = get(options, 'major', 32)
      hector.minor = get(options, 'minor', sqrt3*hector.major/2)
      hector.point_height = get(options, 'point_height', hector.major/4)

      hector.calculate_constraint_factor()
      hector.calculate_inverses()
      hector.calculate_scaling_factors()
   end

   function hector.calculate_constraint_factor()
      -- The constraint factor is used when calculating the 'hex_len'
      -- to determine the correct value of the z-coordinate. It can be
      -- either 1 or -1 depending on the user's choice of x- and
      -- y-axis.
      if hex_to_standard(vector(1, 1)):taxi_len() == 4 then
         hector.constraint_factor = -1
      else
         hector.constraint_factor = 1
      end
   end

   function hector.calculate_inverses()
      -- Not technically inverses but these vectors are used when
      -- converting standard coordinates (3D cube) to user coordinates
      -- (2D hex)
      local x_inv = hector.y_axis:abs()
      x_inv = x_inv / (x_inv..hector.x_axis)
      hector.x_inv = x_inv

      local y_inv = hector.x_axis:abs()
      y_inv = y_inv / (y_inv..hector.y_axis)
      hector.y_inv = y_inv
   end

   function hector.calculate_scaling_factors()
      -- These scaling factors are used when converting to/from screen
      -- coordinates.
      hector.scale_x = (hector.major - hector.point_height ) / 0.866 / sqrt2
      hector.scale_y = hector.minor / sqrt2
   end

   -- Rotation transformations.
   -- [FLAT] and [POINTY] are for converting to/from screen coordinates.
   -- 'sixty' rotates a hex in standard coordinates by sixty degrees.
   hector.rotors = {
      [FLAT]=rotor(0.88047624, 0.11591690, 0.27984814, 0.36470520),
      [POINTY]=rotor(0.88047624, -0.11591690, 0.36470520, 0.27984814),
      sixty=rotor(0.86602540, -0.28867513, 0.28867513, -0.28867513)
   }
   -- END Hector Initialization ------------------------------------------------
end

return hector

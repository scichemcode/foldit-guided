// The Foldit recipe. Kept as a module so the page needs no network request to copy it.
export default String.raw`
-- Foldit move energetics: run once after each move (Cookbook > Load local recipe).
-- Read-only with respect to the protein: no atom is moved. It measures atom-to-anchor
-- distances with temporary disabled bands, deletes them, and prints the result.

local HEAVY = {a=5,r=11,n=8,d=8,c=6,q=9,e=9,g=4,h=10,i=8,l=8,k=9,m=8,f=11,p=7,s=6,t=7,w=14,y=12,v=7}
local MAX_ANCHORS = 10
local baseBands = nil

local function esc(v)
  return (string.gsub(tostring(v or ""), "[^%w%-_]", function(c)
    return string.format("%%%02X", string.byte(c))
  end))
end

local function cleanup(err)
  if baseBands then
    while band.GetCount() > baseBands do band.Delete(band.GetCount()) end
  end
  if undo then undo.SetUndo(true) end
  if err then print("foldit_energy capture stopped: " .. tostring(err)) end
end

-- Foldit refuses bands within a segment and between neighbours: those print as "-".
local function measure(s1, a1, s2, a2)
  if math.abs(s1 - s2) <= 1 then return "-" end
  local before = band.GetCount()
  local b = band.AddBetweenSegments(s1, s2, a1, a2)
  if not b or b == 0 then return "-" end
  local created = band.GetCount() > before
  if created then band.Disable(b) end
  local d = band.GetLength(b)
  if created then band.Delete(b) end
  return string.format("%.6f", d)
end

local function main()
  local ask = dialog.CreateDialog("Capture for energy analysis")
  ask.about = dialog.AddLabel("Takes a snapshot for the energy page. Does not move the protein.")
  ask.move = dialog.AddTextbox("Move", "")
  ask.hint = dialog.AddLabel("Move is optional: what you just did, e.g. shake")
  ask.ok = dialog.AddButton("Capture", 1)
  ask.cancel = dialog.AddButton("Cancel", 0)
  if dialog.Show(ask) ~= 1 then
    print("Capture cancelled. Nothing was recorded.")
    return
  end

  local protein, skipped = {}, 0
  for i = 1, structure.GetCount() do
    local aa = string.lower(structure.GetAminoAcid(i))
    if HEAVY[aa] and structure.GetSecondaryStructure(i) ~= "M" then
      if structure.GetAtomCount(i) < HEAVY[aa] then
        error("segment " .. i .. " reports fewer atoms than expected")
      end
      protein[#protein + 1] = {seg = i, aa = aa}
    else
      skipped = skipped + 1
    end
  end
  local n = #protein
  if n < 11 then error("at least 11 protein segments are needed") end

  -- Alternate CA, O and side-chain tips so a straight chain still gives non-coplanar anchors.
  local k = math.min(MAX_ANCHORS, math.floor((n + 1) / 2))
  local anchors = {}
  for j = 1, k do
    local p = protein[math.floor((j - 1) * (n - 1) / (k - 1)) + 1]
    local choices = {2, 4, HEAVY[p.aa]}
    anchors[j] = {seg = p.seg, atom = choices[(j - 1) % 3 + 1]}
  end

  baseBands = band.GetCount()
  if undo then undo.SetUndo(false) end

  local nonce = tostring(os.time()) .. "-" .. tostring(os.clock())
  print("FEX1|BEGIN|" .. esc(nonce) .. "|" .. esc(puzzle.GetPuzzleID()) .. "|" .. esc(puzzle.GetName())
    .. "|" .. esc(ui.GetTrackName()) .. "|" .. esc(ask.move.value) .. "|" .. esc(current.GetScore())
    .. "|" .. skipped .. "|" .. esc(structure.GetSymCount and structure.GetSymCount() or 0))
  print("FEX1|DESC|" .. esc(puzzle.GetDescription()))

  for _, part in ipairs(puzzle.GetPuzzleSubscoreNames()) do
    local sum = 0
    for i = 1, structure.GetCount() do
      sum = sum + current.GetSegmentEnergySubscore(i, part)
    end
    print("FEX1|PART|" .. esc(part) .. "|" .. esc(sum))
  end

  for j = 1, k do
    local row = {}
    for m = 1, k do
      row[m] = (m == j) and "0" or measure(anchors[j].seg, anchors[j].atom, anchors[m].seg, anchors[m].atom)
    end
    print("FEX1|ANCHOR|" .. j .. "|" .. table.concat(row, ","))
  end

  local count = 0
  local reported = 0
  for index, p in ipairs(protein) do
    -- One atom past the heavy atoms, because a chain end carries an extra OXT.
    for a = 1, math.min(HEAVY[p.aa] + 1, structure.GetAtomCount(p.seg)) do
      local row = {}
      for m = 1, k do row[m] = measure(p.seg, a, anchors[m].seg, anchors[m].atom) end
      print("FEX1|ATOM|" .. p.seg .. "|" .. p.aa .. "|" .. a .. "|" .. table.concat(row, ","))
      count = count + 1
    end
    local percent = math.floor(index * 100 / n)
    if percent < 100 and percent >= reported + 20 then
      reported = percent
      print("Capturing: " .. percent .. "%")
    end
  end

  print("FEX1|END|" .. count)
  cleanup()
  print("Capture finished. Go back to the web page.")
end

xpcall(main, cleanup)
`.trimStart();

%Returns the squared displacements. The first displacement (=0) is dropped.
function sqDisplacements = squaredDisplacements(positions)

nPos = size(positions)(1);

dxs = positions(:,1) - positions(1,1);
dys = positions(:,2) - positions(1,2);
dzs = positions(:,3) - positions(1,3);

sqDisplacements = dxs(2:end).^2 + dys(2:end).^2 + dzs(2:end).^2;


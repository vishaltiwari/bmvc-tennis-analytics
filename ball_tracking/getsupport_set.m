node1_index = 146;
node2_index = 164;

node1 = node_set(node1_index);
node2 = node_set(node2_index);

ss1 = node1.model.support_set;
ss2 = node2.model.support_set;

track = [];
for i=1:size(ss1,2)
    p = ss1(i);
    track = [track ; p.x p.y];
end

for i=2:size(ss2,2)
    p = ss2(i);
    track = [track ; p.x p.y];
end



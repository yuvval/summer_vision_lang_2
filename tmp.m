
c1_det = trk_history(1,2:end);
c2_det = trk_history(2,2:end);

Nframes = size(trk_history,2) -1;
[z1, z1_proj, v1] = deal(nan(1,Nframes));
[z2, z2_proj, v2] = deal(nan(1,Nframes));

feat_name = 'center_z';
feat_id_center_z = find(ismember(tracker_feats.names, feat_name));

feat_name = 'velocity_abs';
feat_id_center_v_abs = find(ismember(tracker_feats.names, feat_name));

for t=1:Nframes
z1(t) = tracker_feats.values{t}(1, c1_det(t), feat_id_center_z);
z2(t) = tracker_feats.values{t}(1, c2_det(t), feat_id_center_z);

v1(t) = tracker_feats.values{t}(1, c1_det(t), feat_id_center_v_abs);
v2(t) = tracker_feats.values{t}(1, c2_det(t), feat_id_center_v_abs);

end

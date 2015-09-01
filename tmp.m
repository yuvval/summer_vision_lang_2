% 
% c1_det = trk_history(1,2:end);
% c2_det = trk_history(2,2:end);
% 
% Nframes = size(trk_history,2) -1;
% [z1, z1_proj, v1] = deal(nan(1,Nframes));
% [z2, z2_proj, v2] = deal(nan(1,Nframes));
% 
% feat_name = 'center_z';
% feat_id_center_z = find(ismember(tracker_feats.names, feat_name));
% 
% feat_name = 'velocity_abs';
% feat_id_center_v_abs = find(ismember(tracker_feats.names, feat_name));
% 
% for t=1:Nframes
% z1(t) = tracker_feats.values{t}(1, c1_det(t), feat_id_center_z);
% z2(t) = tracker_feats.values{t}(1, c2_det(t), feat_id_center_z);
% 
% v1(t) = tracker_feats.values{t}(1, c1_det(t), feat_id_center_v_abs);
% v2(t) = tracker_feats.values{t}(1, c2_det(t), feat_id_center_v_abs);
% 
% end
depth_im_calib = depth_im;
Nframes = length(depth_im_calib);
for t=1:Nframes
    stationary_part = depth_im{t}(15:70,150:200);
    calib_val = mean(stationary_part(:));
    depth_im_calib{t} = depth_im{t} - calib_val;
end

for t=2:Nframes
%     diff_im = depth_im_calib{t} - depth_im_calib{t-1};
    diff_im = depth_im_calib{t};
    imagesc(diff_im.');
%     caxis([0 300]);
    caxis([0 1100]);
%     caxis([-100 100]);
    colorbar;shg
    disp('press any key for the next diff');
    pause

end
disp('done');
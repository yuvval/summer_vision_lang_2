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

load depth_after2sec_behind.mat

depth_im_calib = depth_im;
Nframes = length(depth_im_calib);
for t=1:Nframes
    stationary_part = depth_im{t}(15:70,150:200);
    calib_val = mean(stationary_part(:));
    depth_im_calib{t} = depth_im{t} - calib_val;
end


for t=1:Nframes
    %     show_im = depth_im_calib{t} - depth_im_calib{t-1};
    show_im = depth_im_calib{t};
    imshow(show_im.');
    caxis([0 300]);
    %     caxis([0 1100]);
    %     caxis([-100 100]);
    colormap(flipud(parula));shg
    colorbar;shg
    save_animated_gif_frame('depth_behind.gif', t==1);
    %     disp('press any key for the next diff');
    
end
return
%% 
vid_fname = 'videos/after2sec_behind.avi'; % Person approaches a chair.
frame_sample_interval = 15;
obj = VideoReader(['voc-dpm/' ppvid.vid_fname]);
video = obj.read();

boxes = ppvid.boxes;

trim_first_seconds = 3;
frame_sample_interval = 15;

t=1;
for k=1 + (trim_first_seconds*30):frame_sample_interval:size(video,4)
    
    im=video(:,:,:,k);
    imshow(im);
    save_animated_gif_frame('samples_behind.gif', t==1);
    
    t = t+1;
end
disp('done');

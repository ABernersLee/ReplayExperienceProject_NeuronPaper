function set_ripple_freq(iident)
load(iident,'HP_RippleRaw','CandSeq')
%%%%% gets ripple freq for each cand seq
rp_freq = NaN(size(CandSeq,1),1);
for ii = 1:size(CandSeq,1)
    hp = HP_RippleRaw(HP_RippleRaw(:,1)>=CandSeq(ii,1) & HP_RippleRaw(:,1)<=CandSeq(ii,2),[1 3]);
    [~,locs]=findpeaks(-hp(:,2));  % recording is the reverse LFP
    rp_freq(ii,1) = (length(locs)-1)./(hp(end,1)-hp(1,1));
end

save(iident,'rp_freq','-append')
% Compute the relative perceptual saliency (PS)
function E = Relative_PS(imgIR, imgTV)

E = PS(imgIR)/PS(imgTV);

end
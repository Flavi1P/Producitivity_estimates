From the float 3902681 I want to make oxygen budget profile by profile.

Start simple:
- fixed integration depth (let's say 30m)
- no air-sea fluxes correction
- make a budget with the previous profile
- make budget for night loss (i.e. between a dusk and night profile)
- make budget for net budget (i.e. byeween two dusk or two down profile)
- express the budget in mmol o2 m-2 d-1

You can compare your output with the output stored in @Data\Processed\O2_float_net_change.csv and   @Data\Processed\O2_float_night_loss.csv but I don't know what has been exactly done to produce those estimates. I just know that they are net change and night loss, respectively. 

Once you're happy with the produced estimate make a plot showing the timeseries of your estimates along with the two already existing estimates. 
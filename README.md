# Coral-Dynamic-Energy-Budgets
Symbiotic Dynamic Energy Budget model code for scleractinian corals associated with doi:10.1016/j.jtbi.2009.03.004 and 10.1016/j.ecolmodel.2011.01.004
This is code to solve steady state values of the coral model in Muller et al. (2009) and used in Edmunds et al. (2011). I commented the code reasonably well and there shouldn’t be major bugs or errors (but the code comes ‘as is’ without warranty).
 
The run file runsteadystatecoral2019.m calls the function steadystatecoral.m  Photoinhibion and photodamage modules of Eynaud et al (2011). The plotting section is messy, but you can change this easily to your preference.
 
The file runsteadystatecoralphoto2019 calls steadystatecoralphoto. This version includes the multiplicative linear ROS/linear Damage model of Klansjcek et al (2016). Damage increases maintenance in the symbiont. This combination hasn’t been published. The plotting section is equally messy.
 
There is also code to solve the dynamics of the system. I used this predominantly to check if the steady state solutions made sense.

REFERENCES
Edmunds, P.J., Putnam, H.M., Nisbet, R.M. and Muller, E.B. (2011) Benchmarks in organism performance and their use in comparative analyses. Oecologia 167:379–390. DOI 10.1007/s00442-011-2004-2
Eynaud, Y., Nisbet, R.M. and Muller, E.B. (2011) Impact of excess and harmful radiation on energy budgets in scleractinian corals. Ecological Modelling. 222: 1315-1322. DOI: 10.1016/j.ecolmodel.2011.01.004. 
Klanjšček, T, Muller, E.B. and Nisbet R.M. (2016) Feedbacks and tipping points in organismal response to oxidative stress. Journal of Theoretical Biology 404: 361-374. DOI: 10.1016/j.jtbi.2016.05.034
Muller, E.B. Kooijman, S.A.L.M., Edmunds, P.J. Doyle, F.J. and Nisbet, R.M. 2009. Dynamic energy budgets in syntrophic symbiotic relationships between heterotrophic hosts and photoautotrophic symbionts.  Journal of Theoretical Biology, 259: 44–57. 


# CSEDS_TVMA

### Action Plan

| Date       | Description |
|------------|-------------|
| Thu, 15 Jan | Completed candidate models and data/methodology |
| Sun, 18 Jan | Completed code replication and draft introduction and literature review |
| Thu, 22 Jan | Full analysis and interpretation of results |
| Sun, 25 Jan | Draft of written report and internal review |
| Tue, 27 Jan | Final submission of written report and presentation slides |
| Wed, 28 Jan | Rehearsal and refinement of presentation narrative |
| Thu, 29 Jan | Final presentation and group member evaluation |


### Task Distribution

| Task | Sofia | Marjolein | Naima | Avril |
|------|-------|-----------|-------|-------|
| Model / Data + Methodology |  | A |  | R |
| Code Replication | A |  | R |  |
| Introduction + Literature Review |  | R | A |  |
| Code Extension | A |  |  | R |
| Results Analysis |  | A | R |  |
| Conclusion + Abstract | R |  |  | A |
| Full Report Review |  | R |  | A |
| Slides Preparation | R |  | A |  |

**Legend:** R = responsible, A = assisting

---

**Questions for next Q&A (15 Jan at 14.00)**
- ?

**Team meeting (11 Jan at 16.00)**
- Marjolein and Naima will review overleaf by 21.00 11 JAN
- if anything to fix and improve, me and Avril will stand by till 21.30 
- Submission Monday morning by Marjolein
- Naima not available every Saturday
- Avril and Sofia limited availability during weekdays, very exhaustively evening…
  
**Practical information**
- Overleaf: sofia.e.jeong@gmail.com with password vu202526
- CSEDS Research Proposal - Online LaTeX Editor Overleaf
- Zoom link: https://vu-live.zoom.us/j/98642950060?pwd=eOWIGfakTaQorMoYZ5mK5NXmRTatMa.1

## Research topic
⚠️ Extending weights with structural break-like components might not yield better results, as there is scientific evidence on economic/financial shifts being gradual rather than suddenly happening. 
-> Alternative: The reference paper utilizes fixed-length time windows to determine the weight changes. However, structural changes in the economy / financial market do not occur at equally spaced points in time. Hence, varying length of time windows might be worth looking into!

---

**Avril’s proposal:** Focus on time-varying model averaging (Sun et al., 2021)
-	Apply to economic / financial setting 
-	Research gap: Existing paper proposes model with smooth model weight changes, but what about sudden regime switches? 
-	Extension: Add 2nd weight component with indicator function

**Papers that have extended our paper that we might be able to use:**
-	Sun, Y., Hong, Y., Wang, S., & Zhang, X. (2023). Penalized time-varying model averaging. Journal of Econometrics, 235(2), 1355-1377.
-	Li, H., Zhang, J., Chen, X., & Hong, Y. (2025). TIME-VARYING COMPLETE SUBSET AVERAGING IN A DATA-RICH ENVIRONMENT. Econometric Theory, 1–56. https://doi.org/10.1017/S0266466625000064
-	Sun, Y., Chen, F., & Gao, J. (2024). Model Averaging for Time-Varying Vector Autoregressions. Available at SSRN 5035249.
-	Chen, Q., Hong, Y., & Li, H. (2024). Time-varying forecast combination for factor-augmented regressions with smooth structural changes. Journal of Econometrics, 240(1), 105693.

couldn't find anything on using (Sun et al., 2021) but changing the dgp to have sudden regime changes. So, we could run Monte Carlo simulations that evaluate the TVMA method under sudden regime changes

## Archive

**Naima’s proposal:** Focus on how the choice of information criteria affects the performance of weighted-average regression models (https://doi.org/10.1146/annurev-statistics-041715-033413)
-	Compute model weights using different information criteria
-	Create weighted-average models based on these weights
-	Compare model performance across criteria or look at forecasting (https://doi.org/10.1016/j.ijforecast.2010.04.006)
-	Research gap: comparison of information criteria
-	
**Q&A 09 Jan** 
-	Is it okay to switch from Jackknife criterion to OOS MSE or R2 if too complex to code?
: If Jackknife or any other advanced method does not fit our timeline, we are free to simplify like  SAIC and MSE. 
-	Replicated with the same data (i.e., excess stock returns from Campbell and Thompson’s (2008) popular dataset)?
: We will look for some other papers that cited our ref paper AND provides codes to reproduce. perhaps more empirical perspectives we could gain. 
It's perfectly fine to replicate their finding AND tune some other setting e.g., param, shocks (like we proposed) 
Possible action point to reach out to the authors of the ref paper if they could share the codes

**Team meeting 09 JAN**
-	Waiting for the authors <- just try it out and if we receive, we add it into the paper later (if time allows)
-	Go for one of the alternative AIC and other methods like MSE <- we try this 
-	Look for some our own way like open source github (may not be fully correct)
- AIC, MSE, out-of-sample R2, and possibly Jackknife
-	Main topic: Assumption 9. The bandwidth h = cT−λ for 0 < λ < 1 and 0 < c < ∞.
-	Our own set up for bandwidth: This paper has only considered a global bandwidth for the TVJMA estimator, which may be severely affected by structural changes. It would be preferable to use a time-varying bandwidth for each time point.
-	Main DGP 1 (Smooth Structural Changes)
-	Literature: Dynamic window (bandwidth), Jackknife estimator, TVMA
-	Research proposal
1)	Model (TVMA)
2)	Criterion: AICs and Jackknife(original paper) if Jackknife is not fesasible then we will extend it to the general methods
3)	Impact of the size of bandwidth (+ reference paper)
-	2 main people for writing the draft of the proposal (intro + literature / methodology + research question)  with the other two people review and improve it <- Sofia, Avril by Sunday noon?
-	The other two people will start data prep, look for benchmark code, DGP, and so on. <- Naima, Marjolein 


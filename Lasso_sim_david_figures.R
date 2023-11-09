library(readr)
library(ggplot2)
df_results <- read_csv("David_designs_n20_k7/david_paper_lasso_sim_sampled_ME.csv")



#Graph for main effect sign recovery

q_me_trimmed= ggplot(data = df_results[(df_results$two_FI%in% c(1,4,7)),], aes(x=log(as.numeric(lambda)), y= as.numeric(perfect_me), color= design, linetype=design))+
  geom_line(size=1.2)+
  facet_grid(vars(ME), vars(two_FI))+
  xlab(latex2exp::TeX("$log(\\lambda)"))+
  ylab(latex2exp::TeX("ME $\\phi_{\\lambda}"))+
  #scale_color_manual(values =c("blue","red", "green"))+
  theme_bw()+
  #scale_x_continuous(limits = c(-2,2))+
  # scale_linetype_manual(values =c("solid", "dashed", "solid"))+
  theme( legend.title=element_blank(), text=element_text(size=22))

ggsave("david_lasso_sim_me_trimmed.pdf",q_me_trimmed)



#Graph for 2FI sign recovery


q_2FI_trimmed= ggplot(data = df_results[(df_results$two_FI%in% c(1,4,7)),], aes(x=log(as.numeric(lambda)), y= as.numeric(perfect_2FI), color= design, linetype=design))+
  geom_line(size=1.2)+
  facet_grid(vars(ME), vars(two_FI))+
  xlab(latex2exp::TeX("$log(\\lambda)"))+
  ylab(latex2exp::TeX("2FI $\\phi_{\\lambda}"))+
  #scale_color_manual(values =c("blue","red", "green"))+
  theme_bw()+
  #scale_x_continuous(limits = c(-2,2))+
  # scale_linetype_manual(values =c("solid", "dashed", "solid"))+
  theme( legend.title=element_blank(), text=element_text(size=22))

ggsave("david_lasso_sim_2FI_trimmed.pdf",q_2FI_trimmed)

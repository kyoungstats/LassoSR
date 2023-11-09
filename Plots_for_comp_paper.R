source('/Users/hkyoung/Desktop/SIC_Paper_Code/SIC_Paper_Function_Library.R')

#df_results <- read_csv("PED_Vars_comparison_n14_p24_k7_SN3.csv")
df_results_PED_Var <- read_csv("PED_Vars_comparison_n14_p24_k7_SN3.csv")

df_n12_p26 <- read_csv("marley_woods_comparison_n12_p26.csv")

df_vars <- read_csv("Vars_comparison_efficiency.csv")



q= ggplot(data = df_results_PED_Var, aes(x=log(as.numeric(lambda)), y= as.numeric(P_joint), color=design, linetype = design))+
  geom_line(size=1.2)+
  xlab(latex2exp::TeX("$log(\\lambda)"))+
  ylab(latex2exp::TeX("$\\Phi_{\\lambda}"))+
  scale_color_manual(values =c("blue","red", "green"))+
  theme_bw()+
  scale_x_continuous(limits = c(-2,2))+
  scale_linetype_manual(values =c("solid", "dashed", "solid"))+
  theme( legend.title=element_blank(), text=element_text(size=20))


ggsave("VarsVsPEDn14p24.pdf", q)


q= ggplot(data = df_n12_p26, aes(x=log(as.numeric(lambda)), y= as.numeric(P_joint), color= design, linetype = design))+
  geom_line()+
  # ggtitle("AICc and Sigma Threshold")+
  xlab(latex2exp::TeX("$log(\\lambda)"))+
  ylab(latex2exp::TeX("$\\Phi_{\\lambda}"))+
  scale_color_manual(values =c("blue","red"))+
  # geom_vline(xintercept=log(1.38),linetype = "dotdash", color = "red")+
  # geom_vline(xintercept=log(1.19),linetype = "dotted", color = "red")+
  # geom_vline(xintercept=log(1.25),linetype = "dotdash", color = "blue")+
  # geom_vline(xintercept=log(1.03),linetype = "dotted", color = "blue")+
  # geom_text(label="Fine Grid Lambda", x=-0.05, y =0.005, color="black", angle=90)+
  # geom_text(label="Coarse Grid Lambda", x=0.35, y =0.005, color="black", angle=90)+
  theme_bw()+
  scale_x_continuous(limits = c(-2,2))+
  scale_linetype_manual(values =c("solid", "dashed"))+
  theme( legend.title=element_blank(), text=element_text(size=20), legend.position = c(0.2,.75))

ggsave("marwoodsn12p26.png", dpi=600, device="png", q)


q= ggplot(data = df_vars, aes(x=log(as.numeric(lambda)), y= as.numeric(P_joint), color= design, linetype = design))+
  geom_line()+
  # ggtitle("AICc and Sigma Threshold")+
  xlab(latex2exp::TeX("$log(\\lambda)"))+
  ylab(latex2exp::TeX("$\\Phi_{\\lambda}"))+
  scale_color_manual(values =c("blue","red","green"))+
  
  theme_bw()+
  scale_x_continuous(limits = c(-2,2))+
  scale_linetype_manual(values =c("solid", "dashed","solid"))+
  theme( legend.title=element_blank(), text=element_text(size=20), legend.position = c(0.2,.75))

ggsave("vars_comp.png", dpi=600, device="png", q)

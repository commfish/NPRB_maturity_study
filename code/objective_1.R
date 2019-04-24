# load ---
source('code/helper.r')

# data ----
read.csv("data/prop_pilot_study.csv", check.names = FALSE) %>%
  filter(radcap >0) %>% #delete samples without measurements 
  mutate(rad1 = anu1,#calculate radial units (ie. focus to each annulus)
         rad2 = anu1 + anu2,
         rad3 = anu1+anu2+anu3,
         rad4 = anu1+anu2+anu3+anu4,
         rad5 = anu1+anu2+anu3+anu4+anu5,
         rad6 = anu1+anu2+anu3+anu4+anu5+anu6,
         rad7 = anu1+anu2+anu3+anu4+anu5+anu6+anu7,
         rad8 = anu1+anu2+anu3+anu4+anu5+anu6+anu7+anu8,
         rad9 = anu1+anu2+anu3+anu4+anu5+anu6+anu7+anu8+anu9,
         prop1 = rad1 / radcap,
         prop2 = rad2 / radcap,
         prop3 = rad3 / radcap,
         prop4 = rad4 / radcap,
         prop5 = rad5 / radcap,
         prop6 = rad6 / radcap,
         prop7 = rad7 / radcap,
         prop8 = rad8 / radcap,
         prop9 = rad9 / radcap,
         aprop1 = asintrans(prop1), #transform data by arctransform
         aprop2 = asintrans(prop2),
         aprop3 = asintrans(prop3),
         aprop4 = asintrans(prop4),
         aprop5 = asintrans(prop5),
         aprop6 = asintrans(prop6),
         aprop7 = asintrans(prop7),
         aprop8 = asintrans(prop8),
         aprop9 = asintrans(prop9)) -> data 

# zone data

data %>%
  group_by(zone) %>% 
  mutate(n = n()) %>% 
  group_by(zone, n) %>% 
  summarise_at(vars(contains('prop')), funs(mean, sd), na.rm = T) %>% 
  mutate( se1 = prop1_sd / sqrt(n()),
          se2 = prop2_sd / sqrt(n()),
          se3 = prop3_sd / sqrt(n()),
          se4 = prop4_sd / sqrt(n()),
          se5 = prop5_sd / sqrt(n()),
          se6 = prop6_sd / sqrt(n()),
          se7 = prop7_sd / sqrt(n()),
          se8 = prop8_sd / sqrt(n()),
          se9 = prop9_sd / sqrt(n())) -> data_region

# analysis
# one sample from each fish/zone combo since there are three samples 
# from each fish/zone combo
# backcalculation methods relative to preferred Zone 1

set.seed(167) 

data %>% 
  dplyr::select(id, fish, agecap, lencap, contains('anu'), 
                contains('rad'), zone) %>% 
  group_by(fish, zone) %>% 
  sample_n(1)  %>%  # sample one scale from each fish zone (5 scales sampled from each fish)
  gather(agei, radi, rad1:rad9) %>%
  arrange(id, agei) %>% 
  dplyr::select(id, fish, agecap, lencap, radcap, agei, radi, zone) %>% 
  mutate(agei = as.numeric(stringr::str_sub(agei, 4, 4))) %>% 
  filter(!is.na(radi), agei<=agecap) %>% 
  group_by(zone) %>% 
  nest() %>% 
  mutate(lm_sl = purrr::map(data, ~ lm(radcap ~ lencap, data = .)),
         lm_ls = purrr::map(data, ~ lm(lencap ~ radcap, data = .))) -> fit1

fit1 %>% 
  unnest(lm_sl %>% 
           map(tidy),
         lm_ls %>% 
           map(tidy)) %>% 
  dplyr::select(zone, term, estimate, term1, estimate1) %>%
  mutate(int_sl = case_when(term == '(Intercept)' ~ estimate),
         slope_sl = case_when(term == 'lencap' ~ estimate),
         int_ls = case_when(term1 == '(Intercept)' ~ estimate1),
         slope_ls = case_when(term1 == 'radcap' ~ estimate1)) %>% 
  dplyr::select(-term, -estimate, -term1, -estimate1) %>% 
  group_by(zone) %>% 
  summarise_all(funs(sum), na.rm = T) %>% 
  left_join(fit1 %>% 
              unnest(data)) %>% 
  mutate(SPH_len = (-int_ls / slope_ls) + (lencap + int_ls / slope_ls) * 
           (radi / radcap), 
         BPH_len = lencap * ((int_ls + slope_ls * radi) / (int_ls + slope_ls * 
                                                             radcap))) -> output

output %>% 
  group_by(agei, zone) %>% 
  summarise(n = n(),
            mn_radcap = mean(radcap),
            sd_radcap = sd(radcap),
            ss_radcap = sum(radcap^2),
            cv_radcap = sd_radcap / mn_radcap,
            mn_SPH_len = mean(SPH_len),
            sd_SPH_len = sd(SPH_len),
            ss_SPH_len = sum(SPH_len^2),
            cv_SPH = sd_SPH_len / mn_SPH_len,
            error_SPH = qt(0.975, df= n - 1) * sd_SPH_len / sqrt(n),
            upper_SPH = mn_SPH_len + error_SPH,
            lower_SPH = mn_SPH_len - error_SPH,
            mn_BPH_len = mean(BPH_len),
            sd_BPH_len = sd(BPH_len), 
            ss_BPH_len = sum(BPH_len^2),
            cv_BPH = sd_BPH_len / mn_BPH_len,
            error_BPH = qt(0.975,df = n - 1) * sd_BPH_len / sqrt(n),
            upper_BPH = mn_BPH_len + error_BPH,
            lower_BPH = mn_BPH_len - error_BPH) %T>% 
  write.csv("output/table_summary_obj1.csv")

# Scale proportion hypothesis

output %>% 
  dplyr::select(fish, agei, zone, SPH_len) %>%
  spread(key = zone, value = SPH_len) %>% 
  mutate(A2= abs((A1 - A2) / A1) * 100,
         A3= abs((A1 - A3) / A1) * 100,
         A4= abs((A1 - A4) / A1) * 100,
         A6= abs((A1 - A6) / A1) * 100) %>% 
  dplyr::select(-A1) %>% 
  gather(variable, value, -agei, -fish) %>% 
  group_by(agei, variable) %>% 
  summarise(mean_SPH = mean(value, na.rm=T),
            sd_SPH = sd(value, na.rm=T),
            n_SPH = n(),
            se_SPH = sd(value, na.rm=T) / sqrt(n())) %>%
  mutate(Age = factor(agei),
         Zone = factor(variable)) -> data_wide_sph 

output %>% 
  dplyr::select(fish, agei, zone, BPH_len) %>% 
  spread(zone, BPH_len) %>% 
  mutate(A2= abs((A1 - A2) / A1) * 100,
         A3= abs((A1 - A3) / A1) * 100,
         A4= abs((A1 - A4) / A1) * 100,
         A6= abs((A1 - A6) / A1) * 100) %>% 
  dplyr::select(-A1) %>% 
  gather(variable, value, -agei, -fish) %>% 
  group_by(agei, variable) %>% 
  summarise(mean_BPH = mean(value, na.rm=T),
            sd_BPH = sd(value, na.rm=T),
            n_BPH = n(),
            se_BPH = sd(value, na.rm=T) / sqrt(n())) %>%
  mutate(Age = factor(agei),
         Zone = factor(variable))  -> data_wide_bph 

# plots by SPH
# sph
data_wide_sph %>% 
  group_by(Age) %>% 
  summarise(labels = mean(n_SPH)) -> labels
  
data_wide_sph %>% 
  ggplot(aes(Age, mean_SPH)) +
  geom_bar(aes(fill=Zone), stat="identity", position="dodge", alpha=0.9) +
  scale_fill_grey(start = 0, end = .8, name = "") + 
  theme(legend.position=c(.9,.75)) +
  annotate("text", x = 5, y=11, label="A) Scale Proportional Hypothesis", 
           family="Times") +
  geom_text(data = labels, aes(Age, y = -0.4, label=labels, group=Age), 
            size = 2.5) +
  xlab("Age") +
  ylab("Mean % Difference") -> figa

# bph
data_wide_bph %>% 
  group_by(Age) %>% 
  summarise(labels = mean(n_BPH)) -> labels

data_wide_bph %>% 
  ggplot(aes(Age, mean_BPH)) +
  geom_bar(aes(fill=Zone), stat="identity", position="dodge", alpha=0.9) +
  scale_fill_grey(start = 0, end = .8, guide = F) + 
  theme(legend.position=c(.9,.75)) +
  annotate("text", x = 5, y=30, 
           label="B) Body Proportional Hypothesis", family="Times") +
  geom_text(data = labels, aes(Age, y = -1.2, label=labels, group=Age),
            size = 2.5) +
  xlab("Age") +
  ylab("Mean % Difference") -> figb

x <- grid.arrange(figa, figb, ncol=1)
x
ggsave('figs/length_diff.png', plot = x, height = 8.5, 
         width=6.5, units='in')


# linear model figs ----

output %>% 
  mutate(fit = SPH_len) %>% 
  ggplot(aes(x = lencap, y = radcap)) +
  geom_point(color ="lightgray") +
  geom_line(aes(y = lencap * slope_sl + int_sl)) +
  facet_wrap(~zone) +
  expand_limits(x = c(150, 220)) +
  labs(y = "Scale Radius (mm)", x =  "Capture Length (mm)")

ggsave('figs/SPH_regression.png', height = 6.5, width = 6.5, units = "in")

output %>% 
  mutate(fit = BPH_len) %>% 
  ggplot(aes(x = radcap, y = lencap)) +
  geom_point(color ="lightgray") +
  geom_line(aes(y = radcap * slope_ls + int_ls)) +
  facet_wrap(~zone) +
  labs(x = "Scale Radius (mm)", y =  "Capture Length (mm)") 

  ggsave('figs/BPH_regression.png', height = 6.5, width = 6.5, units = "in")

 
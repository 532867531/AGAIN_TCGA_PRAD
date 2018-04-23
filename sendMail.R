# MyEmail.R
library(mailR)
sender <- "@yeah.net"
recipients <- c("@qq.com")
aletter=send.mail(from = sender,
          to = recipients,
          subject = paste("Program Done.",Sys.timezone(),Sys.time(),sep = "_"),
          body = "My program is finished.",
          smtp = list(host.name = "smtp.yeah.net", port = 25,
                      user.name = "@yeah.net",
                      passwd = "", ssl = FALSE),
          authenticate = TRUE,
          send = FALSE)
r=aletter$send()

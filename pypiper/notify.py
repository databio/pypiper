""" Pypiper notification mechanisms. """

from email.mime.text import MIMEText
from smtplib import SMTP


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



class Mailer(object):
    """ Emailer, defined by the sender. """


    def __init__(self, src):
        """
        Establish the Mailer with the source address.

        :param str src: email address of the sender
        """
        self.src = src


    def send(self, sbj, txt, dest=None):
        """
        Send--with the given subject--the provided message to destination(s).
        If no destinations are given, do a send-to-self.

        :param str sbj: message subject
        :param str txt: message content
        :param str | list[str] dest: where to send, optional;
            if unspecified, send-to-self
        """

        # Send-to-self?
        if not dest:
            dest = self.src

        # Determine the destinations.
        if isinstance(dest, str):
            dests = [dest]
            dest_text = dest
        else:
            # Assume we have an iterable of destinations.
            # Allow TypeError if False.
            dests = list(dest)
            dest_text = ", ".join(dest)

        # Build the message.
        msg = MIMEText(txt)
        msg["Subject"] = sbj
        msg["From"] = self.src
        msg["To"] = dest_text

        # Send the message.
        # TODO: determine how/where to change this.
        s = SMTP("localhost")
        s.sendmail(self.src, dests, msg.as_string())

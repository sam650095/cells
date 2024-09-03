from django.db import models
from django.contrib.auth.models import User

class UploadData(models.Model):
    file_name = models.CharField(max_length=255, blank=True, null=True)

    def __str__(self):
        return f"Step1 File: {self.file_name}"
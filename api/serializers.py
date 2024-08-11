from rest_framework import serializers
import os
from django.conf import settings
class FileUploadSerializer(serializers.Serializer):
    file = serializers.FileField()

    def save(self):
        file = self.validated_data['file']
        
        if file.name.endswith('metadata.csv'):
            target_dir = 'metadata'
        elif file.name.endswith('signal_value.csv'):
            target_dir = 'signal_value'
        else:
            raise serializers.ValidationError('Invalid File Name')

        target_path = os.path.join(settings.MEDIA_ROOT,'tempfile', target_dir)
        if not os.path.exists(target_path):
            os.makedirs(target_path)

        new_file_path = os.path.join(target_path, file.name)
        
        with open(new_file_path, 'wb') as f:
            for chunk in file.chunks():
                f.write(chunk)

        return {'file': new_file_path}
